use popgen::{Count, MultiSiteCounts};
use rust_htslib::bcf::record::GenotypeAllele;
use rust_htslib::bcf::{HeaderRecord, Read};
use std::collections::{BTreeSet, HashMap};
use std::path::Path;
use std::sync::Arc;
use popgen::stats::{GlobalPi, GlobalStatistic};

#[derive(Debug)]
struct Job {
    rid: u32,
    start: u64,
    end: Option<u64>,
}

fn process_job(reader: &mut rust_htslib::bcf::IndexedReader, job: Job) -> (u64, Vec<(Vec<Count>, usize)>) {
    reader.fetch(job.rid, job.start, job.end).unwrap();

    let mut site_counts_from_record = Vec::<popgen::Count>::default();
    let mut ret = vec![];

    for result in reader.records() {
        let record = result.unwrap();

        let num_alleles = record.alleles().len();

        if num_alleles > site_counts_from_record.len() {
            let old_len = site_counts_from_record.len();
            site_counts_from_record.resize(num_alleles, 0);
            site_counts_from_record[..old_len].fill(0);
        } else {
            site_counts_from_record.truncate(num_alleles);
            site_counts_from_record.fill(0);
        }

        let gts = record.genotypes().expect("Error reading genotypes");

        // number of sample in the vcf
        let sample_count = usize::try_from(record.sample_count()).unwrap();
        let mut total_alleles = 0;
        for sample_index in 0..sample_count {
            // for each sample
            for gta in gts.get(sample_index).iter() {
                // Filter out missing genotypes
                if let GenotypeAllele::Phased(ind) | GenotypeAllele::Unphased(ind) = gta {
                    site_counts_from_record[*ind as usize] += 1;
                }
                total_alleles += 1;
            }
        }

        ret.push((site_counts_from_record.clone(), total_alleles));
    }

    (job.start, ret)
}

fn main() {
    let path = Arc::from(Path::new(
        "ALL.chr18.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
    ));
    let bcf = rust_htslib::bcf::IndexedReader::from_path(&path).expect("Error opening file.");
    let header = bcf.header().clone();

    let batch_size = 100_000;

    let contig_lengths = header
        .header_records()
        .into_iter()
        .filter_map(|rec| {
            if let HeaderRecord::Contig { key: _, values } = rec {
                let name = values.get("ID")?;
                let length = values.get("length")?.parse::<u64>().ok()?;
                let rid = header.name2rid(name.as_bytes()).ok()?;
                Some((name.clone(), (rid, length)))
            } else {
                None
            }
        })
        .collect::<HashMap<_, _>>();

    let n_workers = std::thread::available_parallelism().unwrap().get();
    let (tx, rx) = crossbeam_channel::unbounded();
    let (res_tx, res_rx) = crossbeam_channel::unbounded();

    let mut start = 0u64;
    let (rid, contig_len) = contig_lengths["18"];
    while start < contig_len {
        let end = (start + batch_size).min(contig_len);
        tx.send(Job {
            rid,
            start,
            end: Some(end),
        })
            .expect("still alive");
        start = end;
    }

    let handles = (0..n_workers)
        .map(|_| {
            let path = path.clone();
            let rx = rx.clone();
            let res_tx = res_tx.clone();
            std::thread::spawn(move || {
                let mut reader = rust_htslib::bcf::IndexedReader::from_path(&path).unwrap();

                for batch in rx {
                    let r = process_job(&mut reader, batch);
                    res_tx.send(r).expect("still alive");
                }
            })
        })
        .collect::<Vec<_>>();

    drop(tx);
    drop(rx);
    drop(res_tx);

    for h in handles {
        h.join().expect("worker panicked");
    }

    let mut chunks = BTreeSet::new();
    for result in res_rx {
        chunks.insert(result);
    }

    let mut done = MultiSiteCounts::default();
    for site in chunks.into_iter().flat_map(|(start, sites)| sites) {
        done.add_site_from_counts(site.0, site.1 as i32).unwrap();
    }

    dbg!(GlobalPi::try_from_iter_sites(done.iter()));
}

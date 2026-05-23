use noodles::core::{Position, Region};
use noodles::vcf::variant::record::samples::keys::key;
use noodles::vcf::variant::record::samples::series::Value;
use noodles::vcf::variant::record::samples::Sample;
use noodles::vcf::variant::record::AlternateBases;
use noodles::vcf::Header;
use popgen::stats::{GlobalPi, GlobalStatistic};
use popgen::{Count, MultiSiteCounts};
use std::collections::{BTreeSet, HashMap};
use std::path::Path;
use std::sync::Arc;

#[derive(Debug)]
struct Job {
    contig: String,
    header: Arc<Header>,
    start: usize,
    end: Option<usize>,
}

fn process_job<R>(
    reader: &mut noodles::vcf::io::IndexedReader<R>,
    job: Job,
) -> (usize, Vec<(Vec<Count>, u64)>)
where
    R: noodles::bgzf::io::Seek + noodles::bgzf::io::BufRead,
{
    let region = Region::new(
        job.contig.clone(),
        Position::try_from(job.start).unwrap()
            ..=Position::try_from(job.end.unwrap_or(usize::MAX)).unwrap(),
    );
    let query = reader.query(&job.header, &region).unwrap();

    let mut site_counts_from_record = Vec::<popgen::Count>::default();
    let mut ret = vec![];

    let ploidy = 2;

    for result in query.records() {
        let record = result.unwrap();

        let num_alleles = 1 + record.alternate_bases().len();
        if num_alleles > site_counts_from_record.len() {
            let old_len = site_counts_from_record.len();
            site_counts_from_record.resize(num_alleles, 0);
            site_counts_from_record[..old_len].fill(0);
        } else {
            site_counts_from_record.truncate(num_alleles);
            site_counts_from_record.fill(0);
        }

        let mut total_alleles = 0;

        for sample in record.samples().iter() {
            total_alleles += ploidy;
            if let Some(Some(value)) = sample
                // get the GT field
                .get(&job.header, key::GENOTYPE)
                .transpose()
                .unwrap()
            {
                match value {
                    Value::Genotype(genotype) => {
                        for entry in genotype.iter() {
                            site_counts_from_record[entry.unwrap().0.unwrap()] += 1;
                        }
                    }
                    other => {
                        // panic because this has basically no reason to happen
                        dbg!(other);
                        panic!("parsed a genotype field and didn't get a genotype enum variant!");
                    }
                }
            };
        }

        ret.push((site_counts_from_record.clone(), total_alleles));
    }

    (job.start, ret)
}

fn main() {
    let path = Arc::<Path>::from(Path::new(
        "ALL.chr18.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
    ));
    let index = noodles::tabix::fs::read(format!("{}.tbi", (*path).display())).unwrap();

    let mut vcf = noodles::vcf::io::indexed_reader::Builder::default()
        .set_index(index)
        .build_from_path(&path)
        .unwrap();
    let header = Arc::new(vcf.read_header().unwrap());

    let batch_size = 100_000;

    let contig_lengths = header
        .contigs()
        .iter()
        .map(|(name, values)| {
            let length = values.length().unwrap();
            (name.clone(), length)
        })
        .collect::<HashMap<_, _>>();

    let n_workers = std::thread::available_parallelism().unwrap().get();
    let (tx, rx) = crossbeam_channel::unbounded();
    let (res_tx, res_rx) = crossbeam_channel::unbounded();

    let mut start = 1;
    let contig_len = contig_lengths["18"];
    while start <= contig_len {
        let end = (start + batch_size).min(contig_len + 1);
        tx.send(Job {
            header: header.clone(),
            contig: "18".to_string(),
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
                let path = Arc::<Path>::from(Path::new(
                    "ALL.chr18.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
                ));
                let index = noodles::tabix::fs::read(format!("{}.tbi", (*path).display())).unwrap();

                let mut vcf = noodles::vcf::io::indexed_reader::Builder::default()
                    .set_index(index)
                    .build_from_path(&path)
                    .unwrap();

                for batch in rx {
                    let r = process_job(&mut vcf, batch);
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

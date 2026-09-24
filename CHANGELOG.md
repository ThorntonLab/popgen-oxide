# Changelog

All notable changes to this project will be documented in this file.

## [0.10.0-alpha.0] - 2026-09-24

### 🚀 Features

- *(F_ST)* Expose some fields on view struct
- MultiSiteCounts.get
- MutliSiteCounts.{len, is_empty}
- Iterator specialization
- Getters on SiteCounts (#13)
- [**breaking**] Errors (#20)
- Take input data from tree sequences (#21)
- Examples (#29)
- *(tskit)* [**breaking**] Add support for arbitrary sample lists (#100)
- Framework for parallelizable statistic computation (#77)
- [**breaking**] {iter->counts}::SiteCounts::try_new (#114)
- Handle multiple sample set input from a tree sequence (#109)
- *(tskit)* Accept site iterator as input argument (#126)
- [**breaking**] FStatistics semantics and getters (#132)
- [**breaking**] Tskit error variants (#138)
- Windowed counts from tree sequence for single sample set (#136)
- [**breaking**] TryReduce trait (#149)
- Fn to return crate version (#158)

### 🐛 Bug Fixes

- Make iteration robust to empty data
- Iter broke
- Duplicate def
- Fix logic error in when to update right coordinate of current tree (#25)
- *(f_st)* Exclude missing data from denominator (#37)
- *(tskit)* Pass correct value for number of samples (#47)
- Fst now enforces allele ID invariance over pops via new counts type (#62)
- Unwrap try_from_tabular in examples (#91)
- *(tskit)* Correctly handle windows after the last site (#141)
- Tajima's d composition (#150)
- *(tskit)* Correctly handle check for polymorphism across multiple sample sets (#164)

### 📚 Documentation

- Add indexed VCF reading example (#53)
- Deny broken links and set up crate for docs.rs (#118)
- Fix typo in Benjamin Peter's last name (#125)
- Empty stat semantics (#145)
- Change old name in example (#165)

### 🚜 Refactor

- Remove MultiSiteCounts.counts_at
- Use new public len
- Give impl Iterator for MultiSiteCountsIter in trait order (#17)
- Gate noodles and vcf (#15)
- *(F_ST)* Misleading binding name (#36)
- Poor binding name in f_st (#51)
- [**breaking**] Record_to_genotypes_adapter doesn't take num_samples (#50)
- [**breaking**] { -> counts}::MultiSiteCounts (#67)
- [**breaking**] Make GlobalStatistic::add_site fallible (#80)
- [**breaking**] Remove TryFrom impls for stat types (#83)
- [**breaking**] Further downstream fallibility (#84)
- [**breaking**] Remove WhichPopulation (#81)
- [**breaking**] MultiPopulationCounts as linear arrays (#96)
- [**breaking**] Remove re-exports of entire modules (#98)
- Streamline processing of tskit input (#104)
- Only support NodeId for tskit input (#111)
- [**breaking**] Renames (#105)
- *(tskit)* Remove some internal code duplication (#115)
- *(tskit)* Store and manipulate site iterator directly (#116)
- [**breaking**] Rename diversity statsitic (#128)
- FStatistics no longer relies on hashmap (#135)
- [**breaking**] Rename core types (#143)
- [**breaking**] Watterson{->s}Theta (#146)
- [**breaking**] Rename GlobalStatistic -> UnpolarisedSiteStat (#147)
- [**breaking**] Add_site_from_counts from AlleleCounts instead of counts and total (#151)
- Hide AlleleCounts fields, even within the crate (#152)
- [**breaking**] Total_alleles is a Count (i.e. i64) (#154)
- Diversity in one iteration (#153)
- Use slice::as_chunks (#155)
- [**breaking**] Merge count types (#159)
- [**breaking**] Use "sample set" instead of "population" (#160)
- [**breaking**] Erase iterator types (#161)

### 🧪 Testing

- Clippy lints
- Remove unnecessary path prefix (#18)
- Theta_pi no longer feature gated (#19)
- Add some tests involving ancient samples (#23)
- Add test of inline ancient sample (#26)
- Make tskit-based tests run (#45)
- Make sure tskit example works w/o external data (#58)
- Add module to generate random test data (#59)
- Test data now generates explicit genotype arrays (#68)
- Add naive calculations for diversity (pi) (#63)
- Naive implementation of Watterson's theta (#71)
- Add iterator over allele frequencies of 1 for a given number of alleles (#69)
- Add ::testing (#72)
- Naive F_ST (#30)
- Use proptest to test pi (#78)
- *(tskit)* Naive implementations no longer depend on node "sample" status (#102)
- Naive implementation of f2 (#88)
- Add Python integration test module (#130)
- Add explicit rust-side test of reciprocal fixation b/w sample sets (#163)

### ⚙️ Miscellaneous Tasks

- Don't need allow non_snake_case
- Typo fix
- Dead code in tests
- Import ordering in tests
- Update tskit pinned version (#95)
- Set toolchain version for development (#99)
- Manifest info needed for crates.io (#134)

### 💼 Other

- Store pi_S as floats not GlobalPi
- Docs
- Incrementally precompute terms
- Add_population param is "site"
- Take hash set of populations
- Noodles == 0.92.0
- Rand 0.9.0
- Relax noodles constraint (#103)
- Manually specify features (#121)
- Set MSRV (#120)
- Remove crossbeam-channel (#156)

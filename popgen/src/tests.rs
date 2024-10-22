#[cfg(test)]
mod tests {
    use std::iter::repeat_n;
    use rand::seq::SliceRandom;
    use crate::{AlleleCounts, AlleleID};

    #[test]
    fn load_raw() {
        let mut rng = rand::thread_rng();

        let sites = vec![
            {
                let mut site = repeat_n(0, 8)
                    .chain(repeat_n(1, 7))
                    .map(AlleleID::from)
                    .map(Option::from)
                    .collect::<Vec<_>>();
                site.shuffle(&mut rng);
                site
            },
            {
                let mut site = repeat_n(0, 341)
                    .chain(repeat_n(1, 69))
                    .chain(repeat_n(2, 926))
                    .map(AlleleID::from)
                    .map(Option::from)
                    .collect::<Vec<_>>();
                site.shuffle(&mut rng);
                site
            },
        ];

        let counts = AlleleCounts::from_tabular(sites);
        assert_eq!(counts.counts, vec![8, 7, 341, 69, 926]);
        assert_eq!(counts.count_starts, vec![0, 2]);
    }
}

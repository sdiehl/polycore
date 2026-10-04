use polycore::modp_echelon;
use polycore::sample::Rng;
use polycore::{Echelon, Fp, Lead, Modular};

#[test]
fn agrees_with_generic_elimination() {
    for p in [
        2,
        101,
        4_611_686_018_427_387_847,
        18_446_744_073_709_551_557,
    ] {
        for lead in [Lead::Low, Lead::High] {
            let mut rng = Rng::new(19);
            let mut raw = modp_echelon::Echelon::recording(lead, p);
            let mut generic = Echelon::recording(lead);
            for _ in 0..28 {
                // Deliberately include duplicates, zeros, and unordered columns.
                let row: Vec<_> = (0..8)
                    .map(|_| ((rng.next_u64() % 32) as usize, rng.next_u64() % p))
                    .collect();
                let fp: Vec<_> = row.iter().map(|&(j, v)| (j, Fp::new(v, p))).collect();
                assert_eq!(raw.insert(&row), generic.insert(&fp));
            }
            for i in 0..35 {
                let (got, used) = raw.solve(i);
                let (want, dependencies) = generic.solve(i);
                assert_eq!(
                    got,
                    want.iter()
                        .map(|&(j, v)| (j, v.residue_mod(p)))
                        .collect::<Vec<_>>()
                );
                assert_eq!(used, dependencies);
                assert_eq!(raw.uses(i), generic.uses(i));
            }
            let got = raw.nullspace(35);
            let want: Vec<Vec<_>> = generic
                .nullspace(35)
                .iter()
                .map(|r| r.iter().map(|&(j, v)| (j, v.residue_mod(p))).collect())
                .collect();
            assert_eq!(got, want);
        }
    }
}

#[test]
fn rectangular_solve() {
    use polycore::dense::SolveError;
    assert_eq!(
        modp_echelon::solve(&[vec![1, 2], vec![3, 4], vec![4, 6]], &[5, 11, 16], 101),
        Ok(vec![1, 2])
    );
    assert_eq!(
        modp_echelon::solve(&[vec![1, 2], vec![2, 4]], &[5, 11], 101),
        Err(SolveError::Inconsistent)
    );
    assert_eq!(
        modp_echelon::solve(&[vec![1, 2]], &[5], 101),
        Err(SolveError::Deficient)
    );
}

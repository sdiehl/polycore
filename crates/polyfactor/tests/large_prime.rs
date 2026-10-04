use num_traits::{One, Zero};
use polycore::{Fp, Uni};
use polyfactor::factor_mod;

fn check(f: &Uni<Fp>, expected: &[(Uni<Fp>, u32)]) {
    let (lc, got) = factor_mod(f);
    assert_eq!(got.len(), expected.len());
    for part in expected {
        assert!(got.contains(part), "missing {part:?} in {got:?}");
    }
    let product = got
        .iter()
        .fold(Uni::constant(lc), |a, (b, m)| &a * &b.pow(*m));
    assert_eq!(&product, f);
    assert_eq!(factor_mod(f), (lc, got));
}

#[test]
fn large_prime_mixed_degrees_and_multiplicities() {
    for p in [
        101,
        2_305_843_009_213_693_951,
        9_223_372_036_854_775_783,
        18_446_744_073_709_551_557,
    ] {
        let line = |root| Uni::new(vec![-Fp::new(root, p), Fp::new(1, p)]);
        // Find a nonsquare so the quadratic is irreducible.
        let a = (2..p)
            .find(|&a| polycore::modp::pow(a, (p - 1) / 2, p) == p - 1)
            .unwrap();
        let quad = Uni::new(vec![-Fp::new(a, p), Fp::zero(), Fp::one()]);
        let parts = vec![(line(p / 2), 3), (line(p - 2), 1), (line(0), 2), (quad, 2)];
        let f = parts
            .iter()
            .fold(Uni::constant(Fp::new(7, p)), |acc, (g, m)| {
                &acc * &g.pow(*m)
            });
        check(&f, &parts);
    }
}

#[test]
fn small_characteristic_and_constants() {
    for p in [2, 3, 5, 7] {
        let parts: Vec<_> = (0..p)
            .map(|i| {
                (
                    Uni::new(vec![-Fp::new(i, p), Fp::new(1, p)]),
                    u32::try_from(p).unwrap(),
                )
            })
            .collect();
        let f = parts
            .iter()
            .fold(Uni::constant(Fp::new(1, p)), |a, (g, m)| &a * &g.pow(*m));
        check(&f, &parts);
        check(&Uni::constant(Fp::new(1, p)), &[]);
    }
    check(&Uni::zero(), &[]);
}

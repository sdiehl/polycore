//! Berlekamp's algorithm over `GF(p)`: the Frobenius map `v -> v^p` is linear on `GF(p)[x]/(f)`,
//! its fixed points are the constants on each irreducible factor, and `gcd(f, v - s)` for a fixed
//! point `v` and a residue `s` splits `f`.

use num_traits::{One, Zero};
use polycore::dense::nullspace;
use polycore::sample::Rng;
use polycore::{Fp, Uni};

/// The monic irreducible factors of a squarefree monic `f` over `GF(p)`.
pub(crate) fn berlekamp(f: &Uni<Fp>, p: u64) -> Vec<Uni<Fp>> {
    let n = f.deg();
    if n <= 1 {
        return vec![f.clone()];
    }
    let frobenius = Uni::x().powmod(p, f);
    let mut row = Uni::constant(Fp::one());
    let mut m = vec![Vec::new(); n];
    for i in 0..n {
        for (j, col) in m.iter_mut().enumerate() {
            let c = row.0.get(j).copied().unwrap_or_else(Fp::zero);
            col.push(if i == j { c - Fp::one() } else { c });
        }
        row = &(&row * &frobenius) % f;
    }
    let basis = nullspace(&m);
    let mut factors = vec![f.clone()];
    let split = |factors: Vec<Uni<Fp>>, v: &Uni<Fp>| {
        factors
            .into_iter()
            .flat_map(|g| {
                let h = g.gcd(v);
                if h.deg() > 0 && h.deg() < g.deg() {
                    vec![&g / &h, h]
                } else {
                    vec![g]
                }
            })
            .collect::<Vec<_>>()
    };
    if p == 2 {
        // Fixed points take values 0 or 1 on each irreducible component.
        for v in basis.iter().map(|v| Uni::new(v.clone())) {
            factors = split(factors, &v);
        }
    } else {
        // A random fixed point is an independent F_p constant on every component,
        // even when their degrees differ. Its quadratic character partitions them.
        // Each trial needs O(log p) modular polynomial multiplications; there is
        // no enumeration of the field. The seeded Las Vegas loop is reproducible.
        let mut rng = Rng::new(0x6265_726c_656b_616d);
        while factors.len() < basis.len() {
            let mut v = Uni::zero();
            for b in &basis {
                v = &v + &Uni::new(b.clone()).scale(&Fp::new(rng.next_u64() % p, p));
            }
            factors = split(factors, &v);
            let character = &v.powmod((p - 1) / 2, f) - &Uni::constant(Fp::one());
            factors = split(factors, &character);
        }
    }
    factors.sort_by_key(|g| (g.deg(), g.0.iter().map(|c| c.value()).collect::<Vec<_>>()));
    factors
}

/// Squarefree decomposition in characteristic `p`: a vanishing derivative means `f = g(x)^p`.
fn squarefree(f: &Uni<Fp>, p: u64) -> Vec<(Uni<Fp>, u32)> {
    let mut out = Vec::new();
    let mut c = f.gcd(&f.derivative());
    let mut w = f / &c;
    let mut i = 1;
    while w.deg() > 0 {
        let y = w.gcd(&c);
        let z = &w / &y;
        if z.deg() > 0 {
            out.push((z, i));
        }
        w = y;
        c = &c / &w;
        i += 1;
    }
    if c.deg() > 0 {
        let step = usize::try_from(p).unwrap_or(usize::MAX);
        let root = Uni::new(c.0.iter().step_by(step).copied().collect());
        let e = u32::try_from(p).unwrap_or(u32::MAX);
        out.extend(squarefree(&root, p).into_iter().map(|(g, m)| (g, m * e)));
    }
    out
}

/// The leading coefficient and the monic irreducible factors over `GF(p)` with multiplicities,
/// the modulus read from the coefficients: at least one must be bound to a prime `p`. `f =
/// lc * prod g^e`.
///
/// ```
/// use polycore::{Fp, Uni};
/// use polyfactor::factor_mod;
/// let u = |c: &[i64], p| Uni::new(c.iter().map(|&c| Fp::from_i64(c, p)).collect());
/// // 2x^4 + 2 over GF(17): x^4 = -1 has all four roots since 8 divides 16
/// let (lc, parts) = factor_mod(&u(&[2, 0, 0, 0, 2], 17));
/// assert_eq!(lc, Fp::new(2, 17));
/// let roots: Vec<u64> = parts.iter().map(|(g, _)| (-g.0[0]).value()).collect();
/// assert_eq!(roots, [15, 9, 8, 2]);
///
/// // x^2 + 1 is irreducible over GF(7)
/// let (_, parts) = factor_mod(&u(&[1, 0, 1], 7));
/// assert_eq!(parts, vec![(u(&[1, 0, 1], 7), 1)]);
/// ```
pub fn factor_mod(f: &Uni<Fp>) -> (Fp, Vec<(Uni<Fp>, u32)>) {
    if f.is_zero() {
        return (Fp::zero(), Vec::new());
    }
    let p = f.0.iter().map(|c| c.modulus()).max().unwrap_or(0);
    let parts = squarefree(&f.monic(), p)
        .into_iter()
        .flat_map(|(g, m)| berlekamp(&g, p).into_iter().map(move |h| (h, m)))
        .collect();
    (f.lc(), parts)
}

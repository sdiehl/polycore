//! Factoring polynomials in one variable over `GF(p)`, over Q, and over number fields, on top of
//! [`polycore`].
//!
//! # Quickstart
//!
//! Inputs are [`polycore::Uni`]; read one with a [`polycore::Ring`] and [`Uni::from_poly`].
//!
//! ```
//! use polycore::{Order, Ring, Uni};
//! use polyfactor::factor;
//!
//! let ring = Ring::new(["x"], Order::Lex);
//! let f = Uni::from_poly(&ring.parse("x^5 - x^4 - x + 1")?, 0).unwrap();
//! let show = |g: &Uni<_>| ring.show(&g.to_poly(0, 1, Order::Lex));
//! let parts: Vec<(String, u32)> = factor(&f).iter().map(|(g, e)| (show(g), *e)).collect();
//! assert_eq!(parts, [("x + 1".into(), 1), ("x^2 + 1".into(), 1), ("x - 1".into(), 2)]);
//! # Ok::<(), polycore::ParseError>(())
//! ```
//!
//! # Choosing a function
//!
//! | Coefficients | Use | Algorithm |
//! |---|---|---|
//! | `GF(p)`, [`polycore::Fp`] | [`factor_mod`] | Berlekamp |
//! | Q, `BigRational` | [`factor`], [`is_irreducible`] | Hensel lifting and Zassenhaus |
//! | A number field, [`Alg`] | [`factor_over`] | Trager's norm method |
//!
//! # Conventions
//!
//! - Every factor is monic and irreducible, paired with its multiplicity, grouped by the
//!   squarefree part it came from. Over Q and number fields that order is ascending
//!   multiplicity; over `GF(p)` the parts that are `p`-th powers come last.
//! - [`factor`] and [`factor_over`] drop the leading coefficient; [`factor_mod`] returns it.
//! - Polynomials of degree zero have no factors.
//! - The modulus of [`factor_mod`] and the minimal polynomial of a [`NumberField`] are trusted:
//!   neither primality nor irreducibility is checked.
//!
//! [`Uni::from_poly`]: polycore::Uni::from_poly
#![allow(
    clippy::many_single_char_names,
    clippy::must_use_candidate,
    clippy::missing_panics_doc,
    clippy::redundant_pub_crate
)]

mod field;
mod rational;
mod trager;
mod zp;

pub use field::{Alg, NumberField};
pub use rational::{factor, is_irreducible};
pub use trager::factor_over;
pub use zp::factor_mod;

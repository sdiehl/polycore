//! Fields, prime fields, sparse and dense polynomials, and sparse echelon forms: the shared
//! foundation beneath Groebner bases, sparse interpolation and linear system reduction.
//!
//! # Quickstart
//!
//! Read polynomials over Q with a [`Ring`], compute, map the coefficients into a prime field and
//! evaluate there.
//!
//! ```
//! use polycore::{Fp, Order, Ring, combination};
//!
//! let ring = Ring::new(["x", "y"], Order::Lex);
//! let f = ring.parse("x^2*y + x*y^2 + y^2")?;
//! let gens = ring.parse_many("x*y - 1, y^2 - 1")?;
//!
//! let (q, r) = f.divide(&gens);
//! assert_eq!(&combination(&q, &gens) + &r, f);
//! assert_eq!(ring.show(&r), "x + y + 1");
//!
//! let p = 101;
//! let fp = f.map(|c| Fp::from_rational(c, p).unwrap());
//! assert_eq!(fp.eval(&[Fp::new(2, p), Fp::new(3, p)]), Fp::new(39, p));
//! # Ok::<(), polycore::ParseError>(())
//! ```
//!
//! # Choosing a type
//!
//! | Need | Use |
//! |---|---|
//! | Polynomials in several variables, sparse | [`Poly<F>`] with a [`Monomial`] [`Order`] |
//! | Polynomials in one variable, dense | [`Uni<F>`] |
//! | Raw `u64` coefficient slices over one prime, fast | [`fast::PrimeField`] |
//! | Exact rationals | `num_rational::BigRational` (a [`Field`] as is) |
//! | `GF(p)` with `p` chosen at runtime | [`Fp`] |
//! | `GF(P)` with `P` a compile-time constant | [`Gf<P>`] |
//! | Rational functions in one variable | [`RatFunc<F>`] |
//! | Incremental sparse linear algebra | [`Echelon<F>`], or [`modp_echelon`] on bare residues |
//! | Rationals from images modulo many primes | [`crt`] |
//! | Interpolation, recurrences | [`interp::Newton`], [`interp::Thiele`], [`interp::Massey`] |
//!
//! Over a prime field prefer [`Uni::mul_fast`], [`Uni::divrem_fast`] and [`Uni::gcd_fast`] (NTT,
//! Newton division, half-gcd) to `*`, [`Uni::divrem`] and [`Uni::gcd`]: they fall back to
//! schoolbook arithmetic below 64 coefficients, so they are safe to call at any size.
//!
//! # Conventions
//!
//! - [`Field`] is a blanket trait over the operators, so generic code is written with `+ - * /`
//!   on owned values and on `&Poly` / `&Uni` references.
//! - [`Uni`] stores coefficients by *ascending* power, trailing zeros trimmed: `Uni(vec![c0, c1,
//!   c2])` is `c0 + c1 x + c2 x^2`, and the zero polynomial is `Uni(vec![])`.
//! - [`Poly`] stores nonzero terms in *descending* [`Order`], so [`Poly::lt`] is `terms[0]`.
//!   Variables are indices `0..nvars`; a [`Ring`] only attaches names to them.
//! - Operands of a binary operation must share `nvars` and `order`; this is not checked.
//! - [`Uni::deg`] of zero is 0, so test [`Uni::is_zero`] first. [`Uni::gcd`],
//!   [`Uni::bezout`] and the factorizations return monic results.
//! - [`Fp::zero`](num_traits::Zero::zero) and [`Fp::one`](num_traits::One::one) carry no
//!   modulus and adopt the one they meet, so `F::zero()` works in generic code. Mixing two bound
//!   moduli is a bug.
//! - Functions that cannot succeed on some inputs return `Option` (`try_map`, `exact`,
//!   `Fp::from_rational`, `try_inverse`); division by a zero polynomial panics.
//! - Reconstructions in [`crt`] and early termination in [`interp::Massey`] are probabilistic
//!   candidates: verify before trusting them.
//!
//! Not constant time: for scientific computing, never cryptography.
#![allow(
    clippy::missing_errors_doc,
    clippy::missing_panics_doc,
    clippy::many_single_char_names,
    clippy::similar_names,
    clippy::module_name_repetitions,
    clippy::must_use_candidate,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap
)]

pub mod crt;
pub mod dense;
pub mod division;
mod echelon;
pub mod evaluation;
pub mod fast;
mod field;
mod fp;
pub mod interp;
pub mod lehmer;
pub mod modp;
pub mod modp_echelon;
mod monomial;
mod parse;
mod poly;
mod ratfunc;
mod rational;
mod real;
mod residue;
pub mod sample;
mod subgroup;
mod uni;

pub use echelon::{Echelon, Lead, SparseRow};
pub use field::{Field, nat, pow};
pub use fp::{Fp, Gf, Modular, ParseFpError};
pub use modp::Primes;
pub use monomial::{Monomial, Order};
pub use parse::{ParseError, Ring};
pub use poly::{Poly, Term, combination};
pub use ratfunc::RatFunc;
pub use uni::Uni;

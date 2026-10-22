//! `--legacy` against the original Python implementation,
//! `NN_entropy_calc.py`.
//!
//! The expected totals were produced by that script (NumPy 1.26, SciPy 1.13,
//! MDAnalysis 2.8, Python 3.12) on the checked-in fixture, not by this crate.
//! The mode reproduces its arithmetic step for step, so they must agree to the
//! last bit.

use std::path::PathBuf;
use std::process::{Command, Output};

use nn_entropy::bat_library::geometry::{
    bond_angle, legacy_bond_angle, legacy_torsion_angle, torsion_angle,
};
use nn_entropy::bat_library::{CoordinateOptions, Geometry, InternalCoordinates};
use nn_entropy::legacy::legacy_entropy;

/// `python NN_entropy_calc.py test.parm7 -1 test.nc`, all five frames.
const PYTHON_TOTAL_ALL_FRAMES: f64 = -44.640002919814606;
/// `python NN_entropy_calc.py test.parm7 3 test.nc`.
const PYTHON_TOTAL_THREE_FRAMES: f64 = 98.35434987070039;

fn fixture(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(rel)
}

fn run_cli(args: &[&str]) -> Output {
    let top = fixture("tests/fixtures/test.parm7");
    let traj = fixture("tests/fixtures/test.nc");
    Command::new(env!("CARGO_BIN_EXE_nn_entropy"))
        .args([top.to_str().unwrap(), traj.to_str().unwrap()])
        .args(args)
        .output()
        .expect("failed to run CLI")
}

fn cli_entropy(args: &[&str]) -> f64 {
    let output = run_cli(args);
    assert!(output.status.success(), "{args:?} failed");
    String::from_utf8_lossy(&output.stdout)
        .trim()
        .strip_prefix("Total entropy = ")
        .expect("missing entropy output")
        .parse()
        .expect("invalid entropy output")
}

#[test]
fn cli_reproduces_the_python_implementation_bit_for_bit() {
    assert_eq!(cli_entropy(&["--legacy"]), PYTHON_TOTAL_ALL_FRAMES);
    assert_eq!(
        cli_entropy(&["--legacy", "--stop", "3"]),
        PYTHON_TOTAL_THREE_FRAMES
    );
}

#[test]
fn library_reproduces_the_python_implementation_bit_for_bit() {
    let options = CoordinateOptions {
        geometry: Geometry::Legacy,
        ..Default::default()
    };
    let mut internal =
        InternalCoordinates::with_options(&fixture("tests/fixtures/test.parm7"), options).unwrap();
    internal
        .calculate_internal_coords(&fixture("tests/fixtures/test.nc"), usize::MAX, false)
        .unwrap();
    let frames = internal.int_coords.len();
    let data: Vec<Vec<f64>> = (0..internal.dim)
        .map(|j| internal.int_coords.iter().map(|frame| frame[j]).collect())
        .collect();
    assert_eq!(
        legacy_entropy(&data, frames).unwrap(),
        PYTHON_TOTAL_ALL_FRAMES
    );
}

#[test]
fn the_default_estimator_differs_from_the_reference() {
    // Periodic torsions change the result; this is the correction the
    // reference mode exists to undo.
    let default = cli_entropy(&[]);
    assert!((default - PYTHON_TOTAL_ALL_FRAMES).abs() > 1e-3);
}

#[test]
fn options_the_python_implementation_lacked_are_rejected() {
    for flags in [
        &["--exact-constant"][..],
        &["--theiler-window", "1"][..],
        &["--mie-order", "3"][..],
        &["--tree-bonds"][..],
        &["--phase-angles"][..],
        &["--tree-bonds", "--jacobian"][..],
    ] {
        let mut args = vec!["--legacy"];
        args.extend_from_slice(flags);
        let output = run_cli(&args);
        assert!(!output.status.success(), "{args:?} should be rejected");
        assert!(
            String::from_utf8_lossy(&output.stderr).contains("cannot be combined with --legacy"),
            "{args:?}"
        );
    }
}

#[test]
fn legacy_torsions_use_pow_not_sqrt() {
    // NumPy evaluates `x ** 0.5` as pow(x, 0.5), which differs from the
    // correctly rounded sqrt(x) by one ulp for about 1 in 1000 inputs. Over
    // many random torsions the legacy and standard formulas must therefore
    // disagree somewhere, and only at the level of rounding. If the compiler
    // turned pow into sqrt they would never disagree.
    let mut state = 0x2545_f491_4f6c_dd1d_u64;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1u64 << 53) as f64 * 4.0 - 2.0
    };
    let mut differing = 0;
    for _ in 0..20_000 {
        let mut point = || [next(), next(), next()];
        let (a1, a2, a3, a4) = (point(), point(), point(), point());
        let standard = torsion_angle(a1, a2, a3, a4);
        let legacy = legacy_torsion_angle(a1, a2, a3, a4);
        assert!((standard - legacy).abs() <= 1e-14, "{standard} vs {legacy}");
        differing += usize::from(standard != legacy);
    }
    assert!(differing > 0, "legacy torsions never differed from sqrt");
}

#[test]
fn legacy_angles_agree_with_the_standard_formula() {
    let a1 = [1.0, 0.2, -0.1];
    let a2 = [0.0, 0.0, 0.0];
    let a3 = [-0.3, 1.1, 0.4];
    assert!((legacy_bond_angle(a1, a2, a3) - bond_angle(a1, a2, a3)).abs() < 1e-12);
}

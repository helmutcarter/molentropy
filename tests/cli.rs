use std::path::PathBuf;
use std::process::Command;

use nn_entropy::bat_library::{BondSet, CoordinateOptions, InternalCoordinates};
use nn_entropy::{CoordinateMetric, EntropyOptions, FiniteSampleConstant, calculate_entropy};

fn fixture(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(rel)
}

fn parse_entropy(stdout: &[u8]) -> f64 {
    String::from_utf8_lossy(stdout)
        .split('=')
        .next_back()
        .expect("missing entropy output")
        .trim()
        .parse()
        .expect("invalid entropy output")
}

#[test]
fn no_periodic_matches_linear_raw_array_estimator() {
    let top = fixture("tests/fixtures/test.parm7");
    let traj = fixture("tests/fixtures/test.nc");
    let output = Command::new(env!("CARGO_BIN_EXE_nn_entropy"))
        .args([
            top.to_str().unwrap(),
            traj.to_str().unwrap(),
            "--no-periodic",
            "--mie-order",
            "2",
        ])
        .output()
        .expect("failed to run CLI");
    assert!(output.status.success());
    let cli_entropy = parse_entropy(&output.stdout);

    let mut internal = InternalCoordinates::new(&top).unwrap();
    internal
        .calculate_internal_coords(&traj, usize::MAX, false)
        .unwrap();
    let frames = internal.int_coords.len();
    let dimensions = internal.int_coords[0].len();
    let mut data = vec![Vec::with_capacity(frames); dimensions];
    for frame in internal.int_coords {
        for (dimension, value) in frame.into_iter().enumerate() {
            data[dimension].push(value);
        }
    }
    let expected = calculate_entropy(&data, frames, 2, &EntropyOptions::default()).unwrap();
    assert!((cli_entropy - expected).abs() < 1e-12);
}

#[test]
fn periodic_remains_the_default() {
    let top = fixture("tests/fixtures/test.parm7");
    let traj = fixture("tests/fixtures/test.nc");
    let periodic = Command::new(env!("CARGO_BIN_EXE_nn_entropy"))
        .args([top.to_str().unwrap(), traj.to_str().unwrap()])
        .output()
        .unwrap();
    let linear = Command::new(env!("CARGO_BIN_EXE_nn_entropy"))
        .args([
            top.to_str().unwrap(),
            traj.to_str().unwrap(),
            "--no-periodic",
        ])
        .output()
        .unwrap();
    assert!(periodic.status.success() && linear.status.success());
    assert!((parse_entropy(&periodic.stdout) - parse_entropy(&linear.stdout)).abs() > 1.0);
}

#[test]
fn exact_constant_flag_selects_exact_finite_sample_term() {
    let top = fixture("tests/fixtures/test.parm7");
    let traj = fixture("tests/fixtures/test.nc");
    let output = Command::new(env!("CARGO_BIN_EXE_nn_entropy"))
        .args([
            top.to_str().unwrap(),
            traj.to_str().unwrap(),
            "--no-periodic",
            "--exact-constant",
        ])
        .output()
        .unwrap();
    assert!(output.status.success());

    let mut internal = InternalCoordinates::new(&top).unwrap();
    internal
        .calculate_internal_coords(&traj, usize::MAX, false)
        .unwrap();
    let frames = internal.int_coords.len();
    let dimensions = internal.int_coords[0].len();
    let mut data = vec![Vec::with_capacity(frames); dimensions];
    for frame in internal.int_coords {
        for (dimension, value) in frame.into_iter().enumerate() {
            data[dimension].push(value);
        }
    }
    let metrics = vec![CoordinateMetric::Linear; dimensions];
    let expected = calculate_entropy(
        &data,
        frames,
        2,
        &EntropyOptions {
            metrics: Some(&metrics),
            constant: FiniteSampleConstant::Exact,
            ..Default::default()
        },
    )
    .unwrap();
    assert!((parse_entropy(&output.stdout) - expected).abs() < 1e-12);
}

fn run_cli(extra: &[&str]) -> std::process::Output {
    let top = fixture("tests/fixtures/test.parm7");
    let traj = fixture("tests/fixtures/test.nc");
    let mut args = vec![
        top.to_str().unwrap().to_string(),
        traj.to_str().unwrap().to_string(),
    ];
    args.extend(extra.iter().map(|arg| arg.to_string()));
    Command::new(env!("CARGO_BIN_EXE_nn_entropy"))
        .args(args)
        .output()
        .expect("failed to run CLI")
}

#[test]
fn theiler_window_flag_reaches_the_estimator() {
    let unwindowed = run_cli(&[]);
    let zero = run_cli(&["--theiler-window", "0"]);
    let one = run_cli(&["--theiler-window", "1"]);
    assert!(unwindowed.status.success() && zero.status.success() && one.status.success());
    assert_eq!(unwindowed.stdout, zero.stdout);

    let top = fixture("tests/fixtures/test.parm7");
    let traj = fixture("tests/fixtures/test.nc");
    let mut internal = InternalCoordinates::new(&top).unwrap();
    internal
        .calculate_internal_coords(&traj, usize::MAX, false)
        .unwrap();
    let metrics = internal.coordinate_metrics();
    let frames = internal.int_coords.len();
    let mut data = vec![Vec::with_capacity(frames); metrics.len()];
    for frame in internal.int_coords {
        for (dimension, value) in frame.into_iter().enumerate() {
            data[dimension].push(value);
        }
    }
    let expected = calculate_entropy(
        &data,
        frames,
        2,
        &EntropyOptions {
            metrics: Some(&metrics),
            theiler_window: 1,
            ..Default::default()
        },
    )
    .unwrap();
    assert!((parse_entropy(&one.stdout) - expected).abs() < 1e-12);
}

#[test]
fn theiler_window_wider_than_the_trajectory_is_rejected() {
    // The fixture has five frames; w = 2 needs six.
    let output = run_cli(&["--theiler-window", "2"]);
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("Theiler window 2"));
}

#[test]
fn coordinate_flags_reach_the_coordinates_and_the_jacobian() {
    let output = run_cli(&[
        "--tree-bonds",
        "--phase-angles",
        "--jacobian",
        "--stride",
        "2",
    ]);
    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);
    let lines: Vec<&str> = stdout.lines().collect();
    assert_eq!(lines.len(), 2, "unexpected output: {stdout}");
    assert!(lines[0].starts_with("Mean ln J = "));
    assert!(lines[1].starts_with("Total entropy = "));
    let value = |line: &str| line.split('=').next_back().unwrap().trim().parse::<f64>();

    let top = fixture("tests/fixtures/test.parm7");
    let traj = fixture("tests/fixtures/test.nc");
    let options = CoordinateOptions {
        bonds: BondSet::Tree,
        phase_angles: true,
        ..Default::default()
    };
    let mut internal = InternalCoordinates::with_options(&top, options).unwrap();
    internal
        .calculate_internal_coords(&traj, usize::MAX, false)
        .unwrap();
    let metrics = internal.coordinate_metrics();
    let log_jacobian: Vec<f64> = internal
        .log_jacobian()
        .unwrap()
        .into_iter()
        .step_by(2)
        .collect();
    let frames: Vec<&Vec<f64>> = internal.int_coords.iter().step_by(2).collect();
    let mut data = vec![Vec::with_capacity(frames.len()); metrics.len()];
    for frame in &frames {
        for (dimension, value) in frame.iter().enumerate() {
            data[dimension].push(*value);
        }
    }
    let expected_entropy = calculate_entropy(
        &data,
        frames.len(),
        2,
        &EntropyOptions {
            metrics: Some(&metrics),
            ..Default::default()
        },
    )
    .unwrap();
    let expected_log_jacobian = log_jacobian.iter().sum::<f64>() / log_jacobian.len() as f64;
    assert!((value(lines[0]).unwrap() - expected_log_jacobian).abs() < 1e-12);
    assert!((value(lines[1]).unwrap() - expected_entropy).abs() < 1e-12);
}

#[test]
fn default_output_is_only_the_entropy() {
    let output = run_cli(&[]);
    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert_eq!(stdout.lines().count(), 1, "unexpected output: {stdout}");
    assert!(stdout.starts_with("Total entropy = "));
}

#[test]
fn jacobian_without_the_tree_bond_set_is_rejected() {
    for args in [
        &["--jacobian"][..],
        &["--jacobian", "--tree-bonds", "--torsions-only"][..],
    ] {
        let output = run_cli(args);
        assert!(!output.status.success(), "{args:?} should fail");
        assert!(
            String::from_utf8_lossy(&output.stderr).contains("--jacobian requires --tree-bonds")
        );
    }
}

#[test]
fn phase_angles_are_the_default() {
    let default = run_cli(&[]);
    let explicit = run_cli(&["--phase-angles"]);
    let raw = run_cli(&["--no-phase-angles"]);
    assert!(default.status.success() && explicit.status.success() && raw.status.success());
    assert_eq!(default.stdout, explicit.stdout);
    assert_ne!(default.stdout, raw.stdout);

    let top = fixture("tests/fixtures/test.parm7");
    let traj = fixture("tests/fixtures/test.nc");
    let options = CoordinateOptions {
        phase_angles: false,
        ..Default::default()
    };
    let mut internal = InternalCoordinates::with_options(&top, options).unwrap();
    internal
        .calculate_internal_coords(&traj, usize::MAX, false)
        .unwrap();
    let metrics = internal.coordinate_metrics();
    let frames = internal.int_coords.len();
    let mut data = vec![Vec::with_capacity(frames); metrics.len()];
    for frame in internal.int_coords {
        for (dimension, value) in frame.into_iter().enumerate() {
            data[dimension].push(value);
        }
    }
    let expected = calculate_entropy(
        &data,
        frames,
        2,
        &EntropyOptions {
            metrics: Some(&metrics),
            ..Default::default()
        },
    )
    .unwrap();
    assert!((parse_entropy(&raw.stdout) - expected).abs() < 1e-12);
}

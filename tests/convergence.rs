//! The convergence check: thinning, extrapolation, a finite-sample bias with
//! a closed form, and the CLI.

mod common;

use std::f64::consts::TAU;
use std::process::Command;

use common::Rng;
use common::amber::{
    CHAIN_ATOMS, CHAIN_BONDS, NetcdfOptions, Topology, chain_frames, concat_frames, scratch_path,
    write_netcdf, write_prmtop,
};
use nn_entropy::bat_library::InternalCoordinates;
use nn_entropy::convergence::{CONVERGENCE_STRIDES, check_convergence, extrapolate, subsample};
use nn_entropy::split::{Coupling, SoftSelection, SplitOptions, estimate_split_entropy};
use nn_entropy::{CoordinateMetric, EntropyOptions, calculate_entropy};

#[test]
fn thinning_keeps_every_kth_of_the_first_frames() {
    let data = vec![
        (0..10).map(f64::from).collect(),
        (10..20).map(f64::from).collect(),
    ];
    assert_eq!(
        subsample(&data, 9, 4),
        vec![vec![0.0, 4.0, 8.0], vec![10.0, 14.0, 18.0]]
    );
    assert_eq!(subsample(&data, 10, 1), data);
}

#[test]
fn extrapolation_is_exact_for_a_one_over_n_bias() {
    let frames = [250, 500, 1000];
    let at = |n: usize| 5.0 + 30.0 / n as f64;
    let value = extrapolate(&frames, &frames.map(at));
    assert!((value - 5.0).abs() < 1e-12, "{value}");
    assert_eq!(extrapolate(&frames, &[3.0, 3.0, 3.0]), 3.0);
    // A bias falling as 1/sqrt(n) is only partly removed.
    let slow = |n: usize| 1.0 / (n as f64).sqrt();
    let value = extrapolate(&frames, &frames.map(slow));
    assert!(value > 0.0 && value < slow(1000), "{value}");
}

#[test]
fn the_estimate_gets_each_thinned_set_and_its_stride() {
    let data = vec![(0..1001).map(f64::from).collect::<Vec<_>>()];
    let mut calls = Vec::new();
    let check = check_convergence(&data, 1001, |thinned, frames, stride| {
        assert_eq!(thinned[0].len(), frames);
        assert_eq!(thinned[0][1], stride as f64);
        calls.push(stride);
        Ok(1.0 / frames as f64)
    })
    .unwrap();
    assert_eq!(calls, CONVERGENCE_STRIDES);
    let frames: Vec<usize> = check.points.iter().map(|p| p.frames).collect();
    assert_eq!(frames, [251, 501, 1001]);
    assert_eq!(check.change, 1.0 / 1001.0 - 1.0 / 501.0);
    assert!(check.extrapolated.abs() < 1e-12);
    assert!(check.steady());
    let unsteady = check_convergence(&data, 1001, |_, frames, _| {
        Ok(if frames == 501 { 1.0 } else { 0.0 })
    })
    .unwrap();
    assert!(!unsteady.steady());

    let err = check_convergence(&data, 1001, |_, _, stride| {
        if stride == 2 {
            Err("too few".to_string())
        } else {
            Ok(0.0)
        }
    })
    .unwrap_err();
    assert_eq!(err, "every 2nd frame: too few");
    assert!(check_convergence(&data, 2000, |_, _, _| Ok(0.0)).is_err());
}

#[test]
fn the_gaussian_coupling_bias_of_independent_sets_extrapolates_away() {
    // Soft torsions and hard coordinates drawn independently: the coupling is
    // 0, but a Gaussian mutual information from n frames is biased by about
    // -pq / 2n for p soft features (three harmonics, two columns each) and q
    // hard coordinates. The check sees it shrink and extrapolates to 0.
    let mut rng = Rng::new(91);
    let (soft, hard, n) = (3, 20, 4000);
    let mut data: Vec<Vec<f64>> = (0..soft).map(|_| rng.uniform_vec(n, 0.0, TAU)).collect();
    data.extend((0..hard).map(|_| rng.normal_vec(n, 0.0, 1.0)));
    let mut metrics = vec![CoordinateMetric::Periodic { period: TAU }; soft];
    metrics.extend(vec![CoordinateMetric::Linear; hard]);
    let soft_indices: Vec<usize> = (0..soft).collect();
    let settings = SplitOptions::new(SoftSelection::Indices(&soft_indices));
    assert!(matches!(settings.coupling, Coupling::Gaussian { .. }));
    let options = EntropyOptions {
        metrics: Some(&metrics),
        ..Default::default()
    };
    let check = check_convergence(&data, n, |thinned, frames, _| {
        estimate_split_entropy(thinned, frames, &settings, &options).map(|r| r.coupling)
    })
    .unwrap();
    let p = (6 * soft) as f64;
    let q = hard as f64;
    for point in &check.points {
        let expected = -p * q / (2.0 * point.frames as f64);
        assert!(
            (point.entropy - expected).abs() < 0.25 * expected.abs(),
            "{} frames: {} vs {expected}",
            point.frames,
            point.entropy
        );
    }
    let full = check.points[2].entropy;
    let extrapolated = check.extrapolated;
    assert!(
        extrapolated.abs() < 0.4 * full.abs(),
        "{extrapolated} vs {full}"
    );
}

#[test]
fn the_cli_reports_each_thinned_estimate_before_the_total() {
    let mut topology = Topology::default();
    topology.push_molecule("BUT", &CHAIN_ATOMS, &CHAIN_BONDS);
    let top = scratch_path("convergence_chain.parm7");
    let traj = scratch_path("convergence_chain.nc");
    write_prmtop(&top, &topology);
    let mut rng = Rng::new(0xc0ff);
    let frames = chain_frames(&mut rng, 400, [0.0, 0.0, 0.0]);
    write_netcdf(&traj, &concat_frames(&[&frames]), NetcdfOptions::default());

    let run = |args: &[&str]| {
        let output = Command::new(env!("CARGO_BIN_EXE_nn_entropy"))
            .args([top.to_str().unwrap(), traj.to_str().unwrap()])
            .args(args)
            .output()
            .expect("failed to run CLI");
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        String::from_utf8_lossy(&output.stdout)
            .lines()
            .map(str::to_string)
            .collect::<Vec<_>>()
    };
    let lines = run(&["--convergence"]);
    let header = lines
        .iter()
        .position(|l| l == "stride\tframes\tentropy")
        .unwrap();
    let rows: Vec<Vec<&str>> = lines[header + 1..header + 4]
        .iter()
        .map(|l| l.split('\t').collect())
        .collect();
    assert_eq!(rows[0][..2], ["4", "100"]);
    assert_eq!(rows[1][..2], ["2", "200"]);
    assert_eq!(rows[2][..2], ["1", "400"]);
    let total = lines
        .last()
        .unwrap()
        .strip_prefix("Total entropy = ")
        .unwrap();
    assert_eq!(rows[2][2], total);
    assert!(lines[header + 4].starts_with("Change from 200 to 400 frames = "));
    assert!(lines[header + 5].starts_with("Extrapolated entropy = "));

    // The thinned estimates are the library's on every 4th frame.
    let mut internal = InternalCoordinates::new(&top).unwrap();
    internal
        .calculate_internal_coords(&traj, usize::MAX, false)
        .unwrap();
    let data = internal.columns(0, 4).unwrap();
    let metrics = internal.coordinate_metrics();
    let options = EntropyOptions {
        metrics: Some(&metrics),
        ..Default::default()
    };
    let quarter = calculate_entropy(&data, data[0].len(), 2, &options).unwrap();
    assert_eq!(rows[0][2].parse::<f64>().unwrap(), quarter);

    // Every estimator takes the flag; split keeps the full run's soft set.
    for estimator in ["qh", "bqh", "split"] {
        let lines = run(&["--estimator", estimator, "--convergence"]);
        assert!(
            lines
                .iter()
                .any(|l| l.starts_with("Change from 200 to 400"))
        );
    }
}

use std::env;
use std::path::Path;

use nn_entropy::bat_library::{BondSet, CoordinateOptions, Geometry, InternalCoordinates};
use nn_entropy::legacy::legacy_entropy;
use nn_entropy::quasi_harmonic::{
    MarginalEstimator, PeriodicCorrelation, QuasiHarmonicOptions,
    estimate_coordinate_quasi_harmonic_entropy, estimate_quasi_harmonic_entropy,
};
use nn_entropy::{
    CoordinateMetric, EntropyOptions, FiniteSampleConstant, calculate_entropy,
    estimate_coordinate_entropy, estimate_coordinate_mutual_information,
};

#[derive(Clone, Copy, PartialEq, Eq)]
enum Estimator {
    NearestNeighbor,
    QuasiHarmonic,
    BoltzmannQuasiHarmonic,
}

/// Prints `message` and exits with status 1.
fn fail(message: impl std::fmt::Display) -> ! {
    eprintln!("{message}");
    std::process::exit(1);
}

/// The value of `result`, or `context` and the error, then exit with status 1.
fn or_fail<T>(result: Result<T, String>, context: &str) -> T {
    result.unwrap_or_else(|err| fail(format!("{context}: {err}")))
}

/// The value after the flag `args[i]`, parsed as `T`.
fn flag_value<T: std::str::FromStr>(args: &[String], i: usize) -> T {
    let flag = &args[i];
    let Some(raw) = args.get(i + 1) else {
        fail(format!("{flag} requires a value"));
    };
    raw.parse()
        .unwrap_or_else(|_| fail(format!("invalid {flag} value: {raw}")))
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 3 {
        fail(format!(
            "Usage: {} <path_to_parm7> <path_to_nc> [--torsions-only] [--no-periodic] [--exact-constant] [--start N] [--stop N] [--stride N] [--mie-order 1|2|3|4] [--theiler-window W] [--tree-bonds] [--no-phase-angles] [--jacobian] [--legacy] [--estimator nn|qh|bqh] [--histogram-bins N] [--periodic-correlation complex|centered] [--per-coordinate] [--mutual-information]",
            args[0]
        ));
    }
    let top_path = Path::new(&args[1]);
    let traj_path = Path::new(&args[2]);

    let mut torsions_only = false;
    // None until a flag sets it: phase angles are on unless --legacy.
    let mut phase_angles: Option<bool> = None;
    let mut tree_bonds = false;
    let mut jacobian = false;
    let mut no_periodic = false;
    let mut exact_constant = false;
    let mut start: Option<usize> = None;
    let mut stop: Option<usize> = None;
    let mut legacy = false;
    let mut estimator = Estimator::NearestNeighbor;
    let mut histogram_bins: Option<usize> = None;
    let mut per_coordinate = false;
    let mut mutual_information = false;
    let mut periodic_correlation: Option<PeriodicCorrelation> = None;
    let mut mie_order = 2;
    let mut stride = 1usize;
    let mut theiler_window = 0usize;
    let mut i = 3;
    while i < args.len() {
        match args[i].as_str() {
            "--torsions-only" => {
                torsions_only = true;
                i += 1;
            }
            "--phase-angles" => {
                phase_angles = Some(true);
                i += 1;
            }
            "--no-phase-angles" => {
                phase_angles = Some(false);
                i += 1;
            }
            "--tree-bonds" => {
                tree_bonds = true;
                i += 1;
            }
            "--jacobian" => {
                jacobian = true;
                i += 1;
            }
            "--no-periodic" => {
                no_periodic = true;
                i += 1;
            }
            "--exact-constant" => {
                exact_constant = true;
                i += 1;
            }
            "--start" => {
                start = Some(flag_value(&args, i));
                i += 2;
            }
            "--stop" => {
                stop = Some(flag_value(&args, i));
                i += 2;
            }
            "--mie-order" => {
                mie_order = flag_value(&args, i);
                i += 2;
            }
            "--stride" => {
                stride = flag_value(&args, i);
                i += 2;
            }
            "--theiler-window" => {
                theiler_window = flag_value(&args, i);
                i += 2;
            }
            "--legacy" => {
                legacy = true;
                i += 1;
            }
            "--estimator" => {
                if i + 1 >= args.len() {
                    fail("--estimator requires a value");
                }
                estimator = match args[i + 1].as_str() {
                    "nn" => Estimator::NearestNeighbor,
                    "qh" => Estimator::QuasiHarmonic,
                    "bqh" => Estimator::BoltzmannQuasiHarmonic,
                    other => {
                        fail(format!("--estimator must be nn, qh or bqh, got {other}"));
                    }
                };
                i += 2;
            }
            "--mutual-information" => {
                mutual_information = true;
                i += 1;
            }
            "--per-coordinate" => {
                per_coordinate = true;
                i += 1;
            }
            "--periodic-correlation" => {
                if i + 1 >= args.len() {
                    fail("--periodic-correlation requires a value");
                }
                periodic_correlation = Some(match args[i + 1].as_str() {
                    "complex" => PeriodicCorrelation::Complex,
                    "centered" => PeriodicCorrelation::Centered,
                    other => {
                        fail(format!(
                            "--periodic-correlation must be complex or centered, got {other}"
                        ));
                    }
                });
                i += 2;
            }
            "--histogram-bins" => {
                histogram_bins = Some(flag_value(&args, i));
                i += 2;
            }
            other => {
                fail(format!("Unknown argument: {other}"));
            }
        }
    }
    if !(1..=4).contains(&mie_order) {
        fail("--mie-order must be 1, 2, 3, or 4");
    }
    if stride == 0 {
        fail("--stride must be at least 1");
    }
    if jacobian && (!tree_bonds || torsions_only) {
        fail("--jacobian requires --tree-bonds and cannot be combined with --torsions-only");
    }

    if estimator != Estimator::NearestNeighbor && mie_order != 2 {
        fail("--mie-order applies only to --estimator nn");
    }
    if histogram_bins.is_some() && estimator != Estimator::BoltzmannQuasiHarmonic {
        fail("--histogram-bins applies only to --estimator bqh");
    }
    if periodic_correlation.is_some() && estimator == Estimator::NearestNeighbor {
        fail("--periodic-correlation applies only to --estimator qh or bqh");
    }
    if per_coordinate && estimator == Estimator::NearestNeighbor && mie_order > 2 {
        fail("--per-coordinate supports --mie-order 1 or 2");
    }
    if mutual_information && estimator != Estimator::NearestNeighbor {
        fail("--mutual-information applies only to --estimator nn");
    }
    if histogram_bins == Some(0) {
        fail("--histogram-bins must be at least 1");
    }
    if legacy {
        // The original Python estimator had none of these; each would make the
        // result differ from it.
        let conflicts = [
            (exact_constant, "--exact-constant"),
            (theiler_window != 0, "--theiler-window"),
            (mie_order != 2, "--mie-order other than 2"),
            (
                estimator != Estimator::NearestNeighbor,
                "--estimator other than nn",
            ),
            (tree_bonds, "--tree-bonds"),
            (phase_angles == Some(true), "--phase-angles"),
            (jacobian, "--jacobian"),
            (per_coordinate, "--per-coordinate"),
            (mutual_information, "--mutual-information"),
        ];
        for (given, flag) in conflicts {
            if given {
                fail(format!("{flag} cannot be combined with --legacy"));
            }
        }
    }

    let coordinate_options = CoordinateOptions {
        bonds: if tree_bonds {
            BondSet::Tree
        } else {
            BondSet::HeavyAtom
        },
        // The original Python implementation used the raw torsions.
        phase_angles: phase_angles.unwrap_or(!legacy),
        geometry: if legacy {
            Geometry::Legacy
        } else {
            Geometry::Standard
        },
    };
    let mut internal = InternalCoordinates::with_options(top_path, coordinate_options)
        .unwrap_or_else(|err| {
            fail(format!(
                "Cannot read topology {}: {err}",
                top_path.display()
            ))
        });
    internal
        .calculate_internal_coords(traj_path, stop.unwrap_or(usize::MAX), torsions_only)
        .unwrap_or_else(|err| {
            fail(format!(
                "Cannot read trajectory {}: {err}",
                traj_path.display()
            ))
        });

    let start = start.unwrap_or(0);
    let one_d_data = internal
        .columns(start, stride)
        .unwrap_or_else(|err| fail(format!("Cannot select frames: {err}")));
    // Look these up once: the tables below describe every coordinate.
    let kinds = internal.coordinate_kinds();
    let atoms = internal.coordinate_atoms();
    let metrics = if no_periodic {
        vec![CoordinateMetric::Linear; internal.dim]
    } else {
        internal.coordinate_metrics()
    };
    let used_frames = one_d_data[0].len();
    let options = EntropyOptions {
        metrics: Some(&metrics),
        constant: if exact_constant {
            FiniteSampleConstant::Exact
        } else {
            FiniteSampleConstant::Asymptotic
        },
        theiler_window,
    };
    let quasi_harmonic = QuasiHarmonicOptions {
        marginals: match histogram_bins {
            Some(bins) => MarginalEstimator::Histogram { bins },
            None => MarginalEstimator::NearestNeighbor,
        },
        periodic: periodic_correlation.unwrap_or_default(),
    };
    const FAILED: &str = "Entropy calculation failed";

    // The total, and for QH and BQH the first-order and correlation terms it
    // is the sum of.
    let (entropy, terms) = match estimator {
        _ if legacy => (
            or_fail(legacy_entropy(&one_d_data, used_frames), FAILED),
            None,
        ),
        Estimator::NearestNeighbor => (
            or_fail(
                calculate_entropy(&one_d_data, used_frames, mie_order, &options),
                FAILED,
            ),
            None,
        ),
        Estimator::QuasiHarmonic | Estimator::BoltzmannQuasiHarmonic => {
            let result = or_fail(
                estimate_quasi_harmonic_entropy(
                    &one_d_data,
                    used_frames,
                    &quasi_harmonic,
                    &options,
                ),
                FAILED,
            );
            if estimator == Estimator::QuasiHarmonic {
                (
                    result.quasi_harmonic,
                    Some((result.first_order_gaussian, result.correlation)),
                )
            } else {
                (
                    result.boltzmann_quasi_harmonic,
                    Some((result.first_order_boltzmann, result.correlation)),
                )
            }
        }
    };

    let describe = |index: usize| {
        let kind = kinds[index].name();
        // Atoms are numbered from 1, as in Amber and cpptraj.
        let atoms = atoms[index]
            .iter()
            .map(|atom| (atom + 1).to_string())
            .collect::<Vec<_>>()
            .join("-");
        format!("{}\t{kind}\t{atoms}", index + 1)
    };
    if per_coordinate {
        const PER_COORDINATE_FAILED: &str = "Per-coordinate entropy calculation failed";
        if estimator == Estimator::NearestNeighbor {
            let values = or_fail(
                estimate_coordinate_entropy(&one_d_data, used_frames, mie_order, &options),
                PER_COORDINATE_FAILED,
            );
            println!("coordinate\ttype\tatoms\tentropy");
            for (index, value) in values.iter().enumerate() {
                println!("{}\t{value}", describe(index));
            }
        } else {
            let shares = or_fail(
                estimate_coordinate_quasi_harmonic_entropy(
                    &one_d_data,
                    used_frames,
                    &quasi_harmonic,
                    &options,
                ),
                PER_COORDINATE_FAILED,
            );
            let (totals, first_order) = if estimator == Estimator::QuasiHarmonic {
                (&shares.quasi_harmonic, &shares.first_order_gaussian)
            } else {
                (
                    &shares.boltzmann_quasi_harmonic,
                    &shares.first_order_boltzmann,
                )
            };
            println!("coordinate\ttype\tatoms\tentropy\tfirst_order\tcorrelation");
            for index in 0..totals.len() {
                println!(
                    "{}\t{}\t{}\t{}",
                    describe(index),
                    totals[index],
                    first_order[index],
                    shares.correlation[index]
                );
            }
        }
    }
    if mutual_information {
        let values = or_fail(
            estimate_coordinate_mutual_information(&one_d_data, used_frames, &options),
            "Mutual information calculation failed",
        );
        println!("coordinate_a\tcoordinate_b\tmutual_information");
        let pairs = (0..internal.dim).flat_map(|i| ((i + 1)..internal.dim).map(move |j| (i, j)));
        for ((i, j), value) in pairs.zip(&values) {
            println!("{}\t{}\t{value}", i + 1, j + 1);
        }
    }
    if let Some((first_order, correlation)) = terms {
        println!("First-order entropy = {first_order}");
        println!("Correlation term = {correlation}");
    }
    if jacobian {
        // Over the same frames the entropy used.
        let log_jacobian = internal
            .log_jacobian()
            .expect("--tree-bonds without --torsions-only was checked above");
        let mean =
            log_jacobian.iter().skip(start).step_by(stride).sum::<f64>() / used_frames as f64;
        println!("Mean ln J = {}", mean);
    }
    // Last, so scripts that read the final "= value" still get the entropy.
    println!("Total entropy = {}", entropy);
}

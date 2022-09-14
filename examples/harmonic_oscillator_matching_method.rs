use sim_quantum::physics::harmonic_potential;
use sim_quantum::physics::matching::MatchingSolver;
use sim_quantum::utils::{gen_range, relative_error};

use std::fs;

use plotters::prelude::*;

fn main() {
    // Solve the time-independent schrodinger equation using the matching method.
    const STEP_SIZE: f64 = 0.1;
    const INITIAL_ENERGY: f64 = 1.45;
    const INITIAL_ENERGY_STEP_SIZE: f64 = 0.1;
    const ENERGY_STEP_SIZE_CUTOFF: f64 = 0.00001;
    const MIN_X: f64 = -5.0;
    const MAX_X: f64 = 5.0;
    const MATCH_X_VAL: f64 = -1.0;
    let mut using_numerov: bool = true;
    let match_idx = ((MATCH_X_VAL - MIN_X) / STEP_SIZE).round() as usize;
    let mut solver = MatchingSolver::new(
        STEP_SIZE,
        INITIAL_ENERGY,
        INITIAL_ENERGY_STEP_SIZE,
        harmonic_potential,
        ENERGY_STEP_SIZE_CUTOFF,
        MIN_X,
        MAX_X,
        match_idx,
        using_numerov,
    );
    solver.solve();

    // Plot the data
    fs::create_dir_all("img").expect("Failed to create image directory");
    let root_area = BitMapBackend::new("img/harmonic_oscillator_matching_method.png", (1280, 720))
        .into_drawing_area();
    root_area.fill(&WHITE).unwrap();
    let (upper, lower) = root_area.split_vertically((50).percent());

    let mut upper_chart = ChartBuilder::on(&upper)
        .set_label_area_size(LabelAreaPosition::Left, 40)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .caption(
            "Harmonic oscillator wavefunction using the matching method",
            ("sans-serif", 40),
        )
        .build_cartesian_2d(solver.x_min..solver.x_max, -1.0..1.0)
        .unwrap();

    upper_chart.configure_mesh()
        .x_desc("x")
        .y_desc("ψ")
        .axis_desc_style(("sans-serif", 20))
        .draw()
        .unwrap();

    upper_chart
        .draw_series(solver.wavefunction_points().iter().map(|point| {
            Circle::new(
                *point,
                2,
                plotters::style::ShapeStyle {
                    color: if point.0 <= MATCH_X_VAL {
                        BLUE.mix(1.0)
                    } else {
                        RED.mix(1.0)
                    },
                    filled: true,
                    stroke_width: 1,
                },
            )
        }))
        .unwrap()
        .label(format!("E = {:.5}", solver.energy))
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));

    upper_chart
        .draw_series(LineSeries::new(
            gen_range(-10.0_f64..=10.0, 0.01)
                .iter()
                .map(|x| (*x, first_excited_state(*x))),
            &BLACK,
        ))
        .unwrap();

    upper_chart
        .configure_series_labels()
        .position(SeriesLabelPosition::LowerRight)
        .label_font(("sans-serif", 20))
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()
        .unwrap();

    let mut lower_chart = ChartBuilder::on(&lower)
        .set_label_area_size(LabelAreaPosition::Left, 40)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .caption("Relative error (%)", ("sans-serif", 40))
        .build_cartesian_2d(solver.x_min..solver.x_max, -10.0..10.0)
        .unwrap();

    lower_chart.configure_mesh()
        .x_desc("x")
        .y_desc("Relative error (%)")
        .axis_desc_style(("sans-serif", 20))
        .draw()
        .unwrap();

    lower_chart
        .draw_series(LineSeries::new(
            solver.wavefunction_points().iter().map(|point| {
                (
                    point.0,
                    relative_error(point.1, first_excited_state(point.0)) * 100.0,
                )
            }),
            &BLACK,
        ))
        .unwrap()
        .label(if using_numerov { "Numerov method" } else { "Second order second difference method"})
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLACK));

    using_numerov = false;
    let mut solver = MatchingSolver::new(
        STEP_SIZE,
        INITIAL_ENERGY,
        INITIAL_ENERGY_STEP_SIZE,
        harmonic_potential,
        ENERGY_STEP_SIZE_CUTOFF,
        MIN_X,
        MAX_X,
        match_idx,
        using_numerov,
    );
    solver.solve();

    lower_chart
        .draw_series(LineSeries::new(
            solver.wavefunction_points().iter().map(|point| {
                (
                    point.0,
                    relative_error(point.1, first_excited_state(point.0)) * 100.0,
                )
            }),
            &GREEN,
        ))
        .unwrap()
        .label(if using_numerov { "Numerov method" } else { "Second order second difference method"})
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &GREEN));

    lower_chart
        .configure_series_labels()
        .position(SeriesLabelPosition::LowerRight)
        .label_font(("sans-serif", 20))
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()
        .unwrap();
}

fn first_excited_state(x: f64) -> f64 {
    -(1.0 / std::f64::consts::PI).powf(0.25) * 2.0_f64.sqrt() * x * (-0.5 * x * x).exp()
}

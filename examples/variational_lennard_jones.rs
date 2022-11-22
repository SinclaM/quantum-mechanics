use sim_quantum::physics::lennard_jones_potential;
use sim_quantum::physics::variational::VariationalSolver;
use sim_quantum::physics::matching::MatchingSolver;

use std::fs;

use plotters::prelude::*;

fn main() {
    const _STEP_SIZE: f64 = 0.001;
    const _INITIAL_ENERGY: f64 = -5.0;
    const _INITIAL_ENERGY_STEP_SIZE: f64 = 0.1;
    const _ENERGY_STEP_SIZE_CUTOFF: f64 = 0.001;
    const _MIN_X: f64 = 0.5;
    const _MAX_X: f64 = 5.0;
    const _MATCH_X_VAL: f64 = 1.4; 
    const _USING_NUMEROV: bool = true;
    const _GUARDING_SCALE_FACTOR: bool = false;
    let match_idx = ((_MATCH_X_VAL - _MIN_X) / _STEP_SIZE).round() as usize;

    let mut matcher = MatchingSolver::new(
        _STEP_SIZE,
        _INITIAL_ENERGY,
        _INITIAL_ENERGY_STEP_SIZE,
        lennard_jones_potential,
        _ENERGY_STEP_SIZE_CUTOFF,
        _MIN_X,
        _MAX_X,
        match_idx,
        _USING_NUMEROV,
        _GUARDING_SCALE_FACTOR,
    );
    matcher.solve();

    // Solve the time-independent schrodinger equation using the variational Monte-Carlo method.
    const STEP_SIZE: f64 = 0.01;
    const MIN_X: f64 = 0.5;
    const MAX_X: f64 = 5.0;

    let mut solver = VariationalSolver::new(
        STEP_SIZE,
        lennard_jones_potential,
        MIN_X,
        MAX_X,
    );
    solver.wavefunction = matcher.wavefunction_points().iter().map(|(_, psi)| -*psi).collect();
    dbg!(&solver.wavefunction);
    //solver.solve();

    // Plot the data
    fs::create_dir_all("img").expect("Failed to create image directory");
    let root_area = BitMapBackend::new("img/variational_lennard_jones.png", (1280, 720))
        .into_drawing_area();
    root_area.fill(&WHITE).unwrap();

    let psi_min = solver.wavefunction_points().into_iter().map(|(_, psi)| psi).reduce(f64::min).unwrap();
    let psi_max = solver.wavefunction_points().into_iter().map(|(_, psi)| psi).reduce(f64::max).unwrap();
    let mut ctx = ChartBuilder::on(&root_area)
        .set_label_area_size(LabelAreaPosition::Left, 40)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .caption(
            "Wavefunction in a Lennard-Jones potential using the variational Monte-Carlo method",
            ("sans-serif", 40),
        )
        .build_cartesian_2d(solver.x_min..solver.x_max, psi_min..psi_max)
        .unwrap();

    ctx.configure_mesh()
        .x_desc("x")
        .y_desc("ψ")
        .axis_desc_style(("sans-serif", 20))
        .draw()
        .unwrap();

    ctx.draw_series(solver.wavefunction_points().iter().map(|point| {
        Circle::new(
            *point,
            2,
            plotters::style::ShapeStyle {
                color: BLUE.mix(1.0),
                filled: true,
                stroke_width: 1,
            },
        )
    }))
    .unwrap()
    .label(format!("E = {:.3}", solver.energy))
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));

    ctx.configure_series_labels()
        .label_font(("sans-serif", 20))
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()
        .unwrap();
}

use std::fs;

use sim_quantum::prelude::*;
use plotters::prelude::*;

fn main() {
    // Solve the time-independent schrodinger equation using the variational Monte-Carlo method.
    let config = VariationalConfig {
        x_min: 0.5,
        x_max: 5.0,
        step_size: 0.01,
        potential: lennard_jones_potential
    };

    let mut solver = VariationalSolver::new(&config);
    solver.solve();

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
        .build_cartesian_2d(solver.config.x_min..solver.config.x_max, psi_min..psi_max)
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
    .label(format!("E = {:.3}", solver.energy()))
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));

    ctx.configure_series_labels()
        .label_font(("sans-serif", 20))
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()
        .unwrap();
}

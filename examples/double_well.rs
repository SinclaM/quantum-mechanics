use std::fs;

use sim_quantum::prelude::*;
use plotters::prelude::*;

fn main() {
    // Solve the time-independent schrodinger equation using the matching method.
    let config = MatchingConfig {
        x_min: -1.3,
        x_max: 1.3,
        x_match: 0.3,
        step_size: 5e-4,
        initial_energy: 21.0,
        initial_energy_step_size: 1.0,
        energy_step_size_cutoff: 0.001,
        potential: double_well_potential,
        using_numerov: true,
        guarding_scale_factor: true
    };
    let mut solver = MatchingSolver::new(&config);
    solver.solve();

    // Plot the data
    fs::create_dir_all("img").expect("Failed to create image directory");
    let root_area = BitMapBackend::new("img/double_well.png", (1280, 720)).into_drawing_area();
    root_area.fill(&WHITE).unwrap();

    let mut ctx = ChartBuilder::on(&root_area)
        .set_label_area_size(LabelAreaPosition::Left, 40)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .caption(
            "Wavefunction in a double well potential using the matching method",
            ("sans-serif", 40),
        )
        .build_cartesian_2d(config.x_min..config.x_max, -2.0..2.0)
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
                color: if point.0 <= config.x_match {
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
    .label(format!("E = {:.3}", solver.energy()))
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));

    ctx.configure_series_labels()
        .label_font(("sans-serif", 20))
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()
        .unwrap();
}

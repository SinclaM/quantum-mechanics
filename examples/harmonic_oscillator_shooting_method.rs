use std::fs;

use sim_quantum::prelude::*;
use plotters::prelude::*;

fn main() {
    // Solve the time-independent schrodinger equation using the shooting method for
    // odd parity wavefunctions.
    let config = ShootingConfig {
        x_max: 7.0,
        step_size: 0.01,
        initial_energy: 0.0,
        intitial_energy_step_size: 0.01,
        energy_step_size_cutoff: 0.000001,
        wavefunction_cutoff: 100.0,
        potential: harmonic_potential,
        parity: Parity::Odd,
    };

    let mut solver = ShootingSolver::new(&config);
    solver.solve();

    // Plot the data
    fs::create_dir_all("img").expect("Failed to create image directory");
    let root_area = BitMapBackend::new("img/harmonic_oscillator_shooting_method.png", (1280, 720))
        .into_drawing_area();
    root_area.fill(&WHITE).unwrap();

    let mut ctx = ChartBuilder::on(&root_area)
        .set_label_area_size(LabelAreaPosition::Left, 40)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .caption(
            "Harmonic oscillator wavefunction using the shooting method",
            ("sans-serif", 40),
        )
        .build_cartesian_2d(-8.0..8.0, -1.5..1.5)
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

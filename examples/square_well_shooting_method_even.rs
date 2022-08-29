use quantum_mechanics::physics::shooting::{Parity, ShootingSolver};
use quantum_mechanics::physics::{box_potential, L};

use std::fs;

use plotters::prelude::*;

fn main() {
    // Solve the time-independent schrodinger equation using the shooting method for
    // even parity wavefunctions.
    let mut solver =
        ShootingSolver::default(10000, L * 1.15 / 10000.0, 0.0, box_potential, Parity::Even);
    solver.solve();

    // Write the output to a data file
    fs::create_dir_all("img").expect("Failed to create image directory");
    let root_area = BitMapBackend::new("img/square_well_shooting_method_even.png", (1280, 720))
        .into_drawing_area();
    root_area.fill(&WHITE).unwrap();

    let mut ctx = ChartBuilder::on(&root_area)
        .set_label_area_size(LabelAreaPosition::Left, 40)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .caption(
            "Particle-in-a-box wavefunction using the shooting method",
            ("sans-serif", 40),
        )
        .build_cartesian_2d(-1.0..1.0, -0.2..1.2)
        .unwrap();

    ctx.configure_mesh().draw().unwrap();

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

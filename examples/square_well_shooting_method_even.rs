use quantum_mechanics::physics::{box_potential, L};
use quantum_mechanics::physics::shooting::{ShootingSolver, Parity};
use std::fs;
use std::process::Command;

fn main() {
    // Solve the time-independent schrodinger equation using the shooting method for
    // even parity wavefunctions.
    let mut solver = ShootingSolver::default(
        10000,
        L * 1.15 / 10000.0,
        0.0,
        box_potential,
        Parity::Even,
    );
    solver.solve();

    // Write the output to a data file
    fs::create_dir_all("data").expect("Failed to create data directory");
    let mut data_file = fs::File::create("data/square_well_shooting_method_even.txt")
        .expect("Failed to create data file");

    solver
        .dump_to_file(&mut data_file)
        .expect("Failed to write to data file");

    // Plot the data using gnuplot.
    fs::create_dir_all("img").expect("Failed to create image directory");
    Command::new("gnuplot")
        .arg("gnuplot/square_well_shooting_method_even.gpi")
        .output()
        .expect("Failed to run gnuplot");
}

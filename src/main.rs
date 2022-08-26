pub mod physics;

use crate::physics::shooting::ShootingSolverEven;
use std::fs;
use std::process::Command;

fn main() {
    // Solve the time-independent schrodinger equation using the shooting method for
    // even parity wavefunctions.
    let mut solver = ShootingSolverEven::default(
        10000,
        1.1 * physics::L / 10000.0,
        0.0,
        physics::box_potential,
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

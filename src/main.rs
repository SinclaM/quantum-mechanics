pub mod physics;

use crate::physics::shooting::ShootingSolverEven;
use std::fs;
use std::io::Write;
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
    let mut data_file = fs::File::create("data/square_well_shooting_method_even.txt")
        .expect("Failed to create data file");

    let mut x: f64 = -(solver.steps as f64) * solver.step_size;
    for val in solver.wavefunction.iter().rev() {
        write!(data_file, "{} {}\n", x, val).expect("Failed to write to data file");
        x += solver.step_size;
    }
    for val in solver.wavefunction {
        write!(data_file, "{} {}\n", x, val).expect("Failed to write to data file");
        x += solver.step_size;
    }

    // Plot the data using gnuplot.
    fs::create_dir_all("img").unwrap();
    Command::new("gnuplot")
        .arg("gnuplot/square_well_shooting_method_even.gpi")
        .output()
        .expect("Failed to run gnuplot");
}

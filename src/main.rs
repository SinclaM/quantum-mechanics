pub mod physics;

use crate::physics::shooting::ShootingSolverEven;
use f64;
use std::fs;
use std::io::prelude::*;

fn main() {
    const N: usize = 100;
    const DX: f64 = 1.0 / (N as f64);
    const INITIAL_ENERGY: f64 = 0.0;
    const DE: f64 = 10.0;
    const PSI_CUTOFF: f64 = 300.0;
    const ENERGY_STEP_CUTOFF: f64 = 0.0001;

    let mut solver = ShootingSolverEven::new(N, DX, INITIAL_ENERGY, DE, PSI_CUTOFF, box_potential, ENERGY_STEP_CUTOFF);
    solver.solve();
    println!("{:?}", solver.wavefunction);
    println!("{}", solver.energy);

    let mut file = fs::OpenOptions::new().append(true).create(true).open("data.txt").expect("file open failure");
    for ele in solver.wavefunction {
        file.write_all((ele.to_string() + "\n").as_bytes()).expect("file write failure");
    }

}

fn box_potential(x: f64) -> f64 {
    if x.abs() < 1.0 {
        0.0
    } else {
        1000.0
    }
}

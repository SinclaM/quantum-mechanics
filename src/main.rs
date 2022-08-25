pub mod physics;

use crate::physics::shooting::ShootingSolverEven;
use std::fs::File;
use std::io::Write;

fn main() {
    const N: usize = 10000;
    const DX: f64 = 1.0 / (N as f64);
    const INITIAL_ENERGY: f64 = 2.0;
    const DE: f64 = 10.0;
    const PSI_CUTOFF: f64 = 300.0;
    const ENERGY_STEP_CUTOFF: f64 = 0.0001;

    let mut solver = ShootingSolverEven::new(N, DX, INITIAL_ENERGY, DE, PSI_CUTOFF, box_potential, ENERGY_STEP_CUTOFF);
    solver.solve();

    let mut file = File::create("data/square_well_shooting_method.txt").unwrap();
    
    let mut x: f64 = - (N as f64) * DX;
    for ele in solver.wavefunction.iter().rev() {
        write!(file, "{} {}\n", x, ele).unwrap();
        x += DX;
    }
    for ele in solver.wavefunction {
        write!(file, "{} {}\n", x, ele).unwrap();
        x += DX;
    }

}

fn box_potential(x: f64) -> f64 {
    if x.abs() < 1.0 {
        0.0
    } else {
        100000.0
    }
}

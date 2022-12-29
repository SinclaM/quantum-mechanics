pub mod solvers;

pub const L: f64 = 0.5;

pub fn box_potential(x: f64) -> f64 {
    if x.abs() < L {
        0.0
    } else {
        100000.0
    }
}

pub fn double_well_potential(x: f64) -> f64 {
    if x.abs() <= 0.1 {
        100.0
    } else if x.abs() < 1.0 {
        0.0
    } else {
        1e5
    }
}

pub fn harmonic_potential(x: f64) -> f64 {
    0.5 * x * x
}

pub fn lennard_jones_potential(x: f64) -> f64 {
    let epsilon = 10.0;
    let sigma = 1.0;
    4.0 * epsilon * ((sigma / x).powf(12.0) - (sigma / x).powf(6.0))
}

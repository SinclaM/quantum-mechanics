pub mod shooting;
pub mod matching;

pub const L: f64 = 0.5;

pub fn box_potential(x: f64) -> f64 {
    if x.abs() < L {
        0.0
    } else {
        100000.0
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

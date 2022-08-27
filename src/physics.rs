pub mod shooting;

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

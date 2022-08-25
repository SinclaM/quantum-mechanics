pub mod shooting;

pub fn box_potential(x: f64) -> f64 {
    if x.abs() < 1.0 {
        0.0
    } else {
        100000.0
    }
}

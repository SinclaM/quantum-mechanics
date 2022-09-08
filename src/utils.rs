pub mod integration;
pub mod root_finding;


pub fn gen_range(range: std::ops::RangeInclusive<f64>, step: f64) -> Vec<f64> {
    let mut x = *range.start();
    let mut ret = Vec::<f64>::new();
    while x <= *range.end() {
        ret.push(x);
        x += step;
    }
    ret
}

pub fn relative_error(observed: f64, theoretical: f64) -> f64 {
    let mut relative_error = (observed - theoretical) / theoretical;

    if relative_error.abs() > 1e4 {
        relative_error = 1e4 * relative_error.signum();
    }

    relative_error
}

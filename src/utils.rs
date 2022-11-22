pub mod integration;
pub mod finite_difference;
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

pub fn count_nodes(f: &[f64]) -> usize {
    let mut nodes = 0;
    for window in f.windows(2) {
        if window[0].signum() != window[1].signum() {
            nodes += 1;
        }
    }
    nodes
}

#[cfg(test)]
mod tests {
    use crate::utils::{gen_range, count_nodes};
    #[test]
    fn test_count_nodes() {
        let x_vals = gen_range(-2.0..=2.0, 0.01);
        let y_vals: Vec<f64> = x_vals.iter().map(|x| x * x - 1.0).collect();
        assert_eq!(count_nodes(&y_vals), 2);
    }
}


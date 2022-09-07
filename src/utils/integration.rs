pub fn trapezoidal(f: &[f64], delta_x: &f64) -> f64 {
    let mut sum = 0.0;
    for window in f.windows(2) {
        sum += window.iter().sum::<f64>();
    }
    0.5 * delta_x * sum
}

#[cfg(test)]
mod tests {
    use crate::utils::integration::trapezoidal;
    use crate::utils::gen_range;

    #[test]
    fn parabola() {
        let step = 0.01;
        let f: Vec<f64> = gen_range(0.0..=1.0, step)
            .iter_mut()
            .map(|x| *x * *x)
            .collect();
        assert!((0.31..0.35).contains(&trapezoidal(&f, &step)));
    }
}

pub fn trapezoidal(f: &[f64], delta_x: &f64) -> f64 {
    let mut sum = 0.0;
    for window in f.windows(2) {
        sum += window.iter().sum::<f64>();
    }
    0.5 * delta_x * sum
}

#[cfg(test)]
mod tests {
    use std::ops::RangeInclusive;
    use crate::utils::integration::trapezoidal;
    fn gen_range(range: RangeInclusive<f64>, step: f64) -> Vec<f64> {
        let mut x = *range.start();
        let mut ret = Vec::<f64>::new();
        while x <= *range.end() {
            ret.push(x);
            x += step;
        }
        ret
    }

    #[test]
    fn parabola() {
        let step = 0.01;
        let f: Vec<f64> = gen_range(0.0..=1.0, step).iter_mut().map(|x| *x * *x).collect();
        assert!((0.31..0.35).contains(&trapezoidal(&f, &step)));
    }

}

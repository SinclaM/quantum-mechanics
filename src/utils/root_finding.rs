pub fn find_root(
    f: impl Fn(f64) -> f64,
    range: std::ops::RangeInclusive<f64>,
    tolerance: f64,
    max_loops: usize,
) -> Option<f64> {
    let mut x0 = *range.start();
    let mut x1 = *range.end();
    for _ in 1..=max_loops {
        let fx0 = f(x0);
        let fx1 = f(x1);

        let x2 = x1 - fx1 * (x1 - x0) / (fx1 - fx0);
        let interval = (x2 - x1).abs() / x1.abs();

        if interval < tolerance {
            return Some(x2);
        }

        x0 = x1;
        x1 = x2;
    }

    None
}

#[cfg(test)]
mod tests {
    use crate::utils::root_finding::find_root;
    #[test]
    fn test_find_root() {
        let root = find_root(
            |x| (x + 1.0).sqrt() * (0.5 * x).cos().powf(3.0),
            1.0..=2.0,
            1e-15,
            500
        );
        assert!(root.is_some());
        assert!((root.unwrap() - 3.14).abs() < 0.01);
    }
}

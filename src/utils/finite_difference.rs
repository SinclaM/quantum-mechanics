pub enum SecondDerivateMethod {
    ForwardDifference,
    CentralDifference,
    BackwardDifference,
}

pub fn second_derivative(
    method: &SecondDerivateMethod,
    f: &[f64],
    i: usize,
    step_size: f64,
) -> f64 {
    match method {
        SecondDerivateMethod::ForwardDifference => {
            (f[i + 2] - 2.0 * f[i + 1] + f[i]) / (step_size * step_size)
        }
        SecondDerivateMethod::CentralDifference => {
            (f[i + 1] - 2.0 * f[i] + f[i - 1]) / (step_size * step_size)
        }
        SecondDerivateMethod::BackwardDifference => {
            (f[i] - 2.0 * f[i - 1] + f[i - 2]) / (step_size * step_size)
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::utils::finite_difference::*;
    use crate::utils::{gen_range, relative_error};
    #[test]
    fn second_derivative_test() {
        let f: Vec<f64> = gen_range(0.0..=1.0, 0.01)
            .iter_mut()
            .map(|x| *x * *x)
            .collect();
        assert!(
            relative_error(
                2.0,
                second_derivative(&SecondDerivateMethod::ForwardDifference, &f, 0, 0.01)
            )
            .abs()
                < 0.01
        );

        assert!(
            relative_error(
                2.0,
                second_derivative(&SecondDerivateMethod::CentralDifference, &f, 50, 0.01)
            )
            .abs()
                < 0.01
        );

        assert!(
            relative_error(
                2.0,
                second_derivative(&SecondDerivateMethod::BackwardDifference, &f, 99, 0.01)
            )
            .abs()
                < 0.01
        );
    }

}

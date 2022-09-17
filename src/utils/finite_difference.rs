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
        },
        SecondDerivateMethod::BackwardDifference => {
            (f[i] - 2.0 * f[i - 1] + f[i - 2]) / (step_size * step_size)
        },
    }
}

use crate::utils::finite_difference;
use crate::utils::finite_difference::SecondDerivateMethod;
use crate::utils::integration;

pub struct VariationalSolver {
    pub steps: usize,
    pub step_size: f64,
    pub last_energy: Option<f64>,
    pub potential: fn(f64) -> f64,
    pub x_min: f64,
    pub x_max: f64,
    pub wavefunction: Vec<f64>,
}

impl VariationalSolver {
    pub fn new(step_size: f64, potential: fn(f64) -> f64, x_min: f64, x_max: f64) -> Self {
        let steps = ((x_max - x_min) / step_size).round() as usize + 1;
        let mut wavefunction = vec![(1.0 / (x_max - x_min)).sqrt(); steps];
        for (i, psi) in wavefunction.iter_mut().enumerate() {
            *psi = (i as f64).powf(2.0) * 4e-7;
        }

        VariationalSolver {
            steps,
            step_size,
            last_energy: None,
            potential,
            x_min,
            x_max,
            wavefunction,
        }
    }

    pub fn x_from_index(&self, i: usize) -> f64 {
        self.x_min + (i as f64) * self.step_size
    }

    fn hamiltonian_on_wavefunction(&self) -> Vec<f64> {
        let mut result = Vec::with_capacity(self.steps);

        for i in 0..=(self.steps - 1) {
            let method = if i == 0 {
                SecondDerivateMethod::ForwardDifference
            } else if i == self.steps - 1 {
                SecondDerivateMethod::BackwardDifference
            } else {
                SecondDerivateMethod::CentralDifference
            };
            result.push(
                -0.5 * finite_difference::second_derivative(
                    &method,
                    &self.wavefunction,
                    i,
                    self.step_size,
                ),
            );
        }
        result
    }

    pub fn current_energy(&self) -> f64 {
        let mut psi_hamil_psi = Vec::with_capacity(self.steps);
        let mut psi_psi = Vec::with_capacity(self.steps);

        for (psi, hamiltonian_on_psi) in self
            .wavefunction
            .iter()
            .zip(self.hamiltonian_on_wavefunction().iter())
        {
            psi_hamil_psi.push(psi * hamiltonian_on_psi);
            psi_psi.push(psi * psi);
        }

        integration::trapezoidal(&psi_hamil_psi, &self.step_size)
            / integration::trapezoidal(&psi_psi, &self.step_size)
    }

    pub fn wavefunction_points(&self) -> Vec<(f64, f64)> {
        let mut pairs: Vec<(f64, f64)> = Vec::with_capacity(self.steps);

        for (i, psi_val) in self.wavefunction.iter().enumerate() {
            pairs.push((self.x_from_index(i), *psi_val));
        }

        pairs
    }
}

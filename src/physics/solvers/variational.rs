use crate::physics::solvers::Solver;
use crate::utils::finite_difference;
use crate::utils::finite_difference::SecondDerivateMethod;
use crate::utils::integration;

use rand::Rng;

#[derive(Clone)]
pub struct VariationalConfig {
    pub x_min: f64,
    pub x_max: f64,
    pub step_size: f64,
    pub potential: fn(f64) -> f64,
}

pub struct VariationalSolver {
    pub config: VariationalConfig,
    steps: usize,
    energy: f64,
    wavefunction: Vec<f64>,
    last_energy: Option<f64>,
}

impl VariationalSolver {
    fn step(&mut self) {
        let mut candidate: Vec<f64> = self.wavefunction.iter().cloned().collect();

        let index = rand::thread_rng().gen_range(0, self.steps);

        let max_delta = 0.01;
        let psi_delta = rand::thread_rng().gen_range(-max_delta, max_delta);
        candidate[index] += psi_delta;

        let candidate_energy = energy_of(
            &candidate,
            self.steps,
            self.config.step_size,
            self.config.potential,
            self.config.x_min,
        );
        if candidate_energy < self.energy {
            self.last_energy = Some(self.energy);

            //println!(
            //"At x = {}: {} -> {}",
            //x_from_index(index, self.x_min, self.step_size),
            //self.wavefunction[index],
            //candidate[index]
            //);
            self.energy = candidate_energy;
            self.wavefunction = candidate;
            println!("{}", self.energy);
        }
    }

    fn normalize(&mut self) {
        let (_, mut f): (Vec<_>, Vec<_>) = self.wavefunction_points().iter().cloned().unzip();
        f = f.iter().map(|val| val * val).collect();
        let integral = integration::trapezoidal(&f, &self.config.step_size);
        self.wavefunction
            .iter_mut()
            .for_each(|val| *val = *val * (1.0 / integral).sqrt());
    }

    pub fn wavefunction_points(&self) -> Vec<(f64, f64)> {
        let mut pairs: Vec<(f64, f64)> = Vec::with_capacity(self.steps);

        for (i, psi_val) in self.wavefunction.iter().enumerate() {
            pairs.push((x_from_index(i, self.config.x_min, self.config.step_size), *psi_val));
        }

        pairs
    }
}

impl Solver for VariationalSolver {
    type CONFIG = VariationalConfig;

    fn new(config: &Self::CONFIG) -> Self {
        let steps = ((config.x_max - config.x_min) / config.step_size).round() as usize + 1;
        let mut wavefunction = vec![(1.0 / (config.x_max - config.x_min)).sqrt(); steps];
        for (i, val) in wavefunction.iter_mut().enumerate() {
            let x = x_from_index(i, config.x_min, config.step_size);
            if x < -1.0 || x > 1.0 {
                *val = 0.0;
            }
        }
        //for (i, psi) in wavefunction.iter_mut().enumerate() {
        //*psi = (i as f64).powf(2.0) * 4e-7;
        //}
        VariationalSolver {
            config: config.clone(),
            steps,
            energy: energy_of(
                &wavefunction,
                steps,
                config.step_size,
                config.potential,
                config.x_min,
            ),
            wavefunction,
            last_energy: None,
        }
    }

    fn solve(&mut self) {
        for i in 1..=100000 {
            if i % 10000 == 0 {
                println!("{i}, {}", self.energy)
            };
            self.step();
        }
        self.normalize();
        //self.energy = energy_of(
        //&self.wavefunction,
        //self.steps,
        //self.step_size,
        //self.potential,
        //self.x_min,
        //);
    }

    fn reset(&mut self) {
        *self = Self::new(&self.config);
    }

    fn energy(&self) -> f64 {
        energy_of(
            &self.wavefunction,
            self.steps,
            self.config.step_size,
            self.config.potential,
            self.config.x_min,
        )
    }

    fn wavefunction_points(&self) -> Vec<(f64, f64)> {
        let mut pairs: Vec<(f64, f64)> = Vec::with_capacity(self.steps);

        for (i, psi_val) in self.wavefunction.iter().enumerate() {
            pairs.push((x_from_index(i, self.config.x_min, self.config.step_size), *psi_val));
        }

        pairs
    }
}

fn x_from_index(i: usize, x_min: f64, step_size: f64) -> f64 {
    x_min + (i as f64) * step_size
}

fn hamiltonian_on_wavefunction(
    wavefunction: &[f64],
    steps: usize,
    step_size: f64,
    potential: fn(f64) -> f64,
    x_min: f64,
) -> Vec<f64> {
    let mut result = Vec::with_capacity(steps);

    for i in 0..=(steps - 1) {
        let method = if i == 0 {
            SecondDerivateMethod::ForwardDifference
        } else if i == steps - 1 {
            SecondDerivateMethod::BackwardDifference
        } else {
            SecondDerivateMethod::CentralDifference
        };
        result.push(
            -0.5 * finite_difference::second_derivative(&method, &wavefunction, i, step_size)
                + potential(x_from_index(i, x_min, step_size)) * wavefunction[i],
        );
    }
    result
}

pub fn energy_of(
    wavefunction: &[f64],
    steps: usize,
    step_size: f64,
    potential: fn(f64) -> f64,
    x_min: f64,
) -> f64 {
    let mut psi_hamil_psi = Vec::with_capacity(steps);
    let mut psi_psi = Vec::with_capacity(steps);

    for (psi, hamiltonian_on_psi) in wavefunction
        .iter()
        .zip(hamiltonian_on_wavefunction(wavefunction, steps, step_size, potential, x_min).iter())
    {
        psi_hamil_psi.push(psi * hamiltonian_on_psi);
        psi_psi.push(psi * psi);
    }

    integration::trapezoidal(&psi_hamil_psi, &step_size)
        / integration::trapezoidal(&psi_psi, &step_size)
}

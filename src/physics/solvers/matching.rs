//! Matching method for solving the time-independent Schrodinger equation
//! in one dimension.

use crate::physics::solvers::Solver;
use crate::utils::integration::trapezoidal;

#[derive(Clone)]
pub struct MatchingConfig {
    pub x_min: f64,
    pub x_max: f64,
    pub x_match: f64,
    pub step_size: f64,
    pub initial_energy: f64,
    pub initial_energy_step_size: f64,
    pub energy_step_size_cutoff: f64,
    pub potential: fn(f64) -> f64,
    pub using_numerov: bool,
    pub guarding_scale_factor: bool,
}

/// A solver that looks for solutions using the matching method.
pub struct MatchingSolver {
    pub config: MatchingConfig,
    steps: usize,
    energy: f64,
    energy_step_size: f64,
    is_left_slope_larger: Option<bool>,
    left_wavefunction: Vec<f64>,
    right_wavefunction: Vec<f64>,
}

impl MatchingSolver {
    fn match_idx(&self) -> usize {
        ((self.config.x_match - self.config.x_min) / self.config.step_size).round() as usize
    }

    /// Returns the x value associated with an index into either the right or left
    /// wavefunction vector.
    fn x_from_index(&self, i: usize, side: &Side) -> f64 {
        match side {
            Side::Left => self.config.x_min + (i as f64) * self.config.step_size,
            Side::Right => self.config.x_max - (i as f64) * self.config.step_size,
        }
    }

    /// Computes a term needed for the Numerov method.
    fn k_sqr(&self, x: f64) -> f64 {
        2.0 * (self.energy - (self.config.potential)(x))
    }

    /// Applies the finite difference approximation to find the value of wavefunction
    /// one position toward the matching point from either the left or right.
    fn next(&self, side: &Side, last_index: usize, psi_last: f64, psi_second_to_last: f64) -> f64 {
        let next: f64;
        if self.config.using_numerov {
            next = (2.0
                * (1.0
                    - (5.0 / 12.0)
                        * self.config.step_size.powf(2.0)
                        * self.k_sqr(self.x_from_index(last_index, &side)))
                * psi_last
                - (1.0
                    + (1.0 / 12.0)
                        * self.config.step_size.powf(2.0)
                        * self.k_sqr(self.x_from_index(last_index - 1, &side)))
                    * psi_second_to_last)
                / (1.0
                    + (1.0 / 12.0)
                        * self.config.step_size.powf(2.0)
                        * self.k_sqr(self.x_from_index(last_index + 1, &side)));
        } else {
            next = 2.0
                * (self.config.step_size
                    * self.config.step_size
                    * ((self.config.potential)(self.x_from_index(last_index, &side))
                        - self.energy)
                    + 1.0)
                * psi_last
                - psi_second_to_last;
        }
        next
    }

    /// Updates the corresponding wavefunction according to the finite difference approximation
    fn step(&mut self, side: &Side) {
        match side {
            Side::Left => {
                let last_index = self.left_wavefunction.len() - 1;
                self.left_wavefunction.push(self.next(
                    side,
                    last_index,
                    self.left_wavefunction[last_index],
                    self.left_wavefunction[last_index - 1],
                ))
            }
            Side::Right => {
                let last_index = self.right_wavefunction.len() - 1;
                self.right_wavefunction.push(self.next(
                    side,
                    last_index,
                    self.right_wavefunction[last_index],
                    self.right_wavefunction[last_index - 1],
                ))
            }
        }
    }

    /// Approximates the wavefunction for the current energy. Stops when it has
    /// computed the requested number of steps, or if the wavefunction begins
    /// diverging.
    fn compute_wavefunction(&mut self) -> Result<(), ()> {
        self.reset_wavefunction();
        for _ in 1..=(self.match_idx() - 1) {
            self.step(&Side::Left);
        }

        for _ in 1..=(self.steps - self.match_idx() - 2) {
            self.step(&Side::Right);
        }

        // Ensure continuity at match_idx
        let scale_factor =
            self.left_wavefunction.last().unwrap() / self.right_wavefunction.last().unwrap();

        if self.config.guarding_scale_factor
            && (scale_factor.abs() > 100.0 || scale_factor.abs().recip() > 100.0)
        {
            Err(())
        } else {
            self.right_wavefunction
                .iter_mut()
                .for_each(|val| *val *= scale_factor);
            Ok(())
        }
    }

    /// Resets the wavefunction vectors.
    fn reset_wavefunction(&mut self) {
        self.left_wavefunction.clear();
        self.right_wavefunction.clear();

        self.left_wavefunction.push(0.0);
        self.left_wavefunction.push(1e-18 * self.config.step_size);

        self.right_wavefunction.push(0.0);
        self.right_wavefunction.push(1e-18 * self.config.step_size);
    }

    /// Returns the slopes of both component wavefunctions (left and right) at the
    /// matching point as a tuple, `(left_slope, right_slope)`.
    pub fn slopes(&self) -> (f64, f64) {
        let left_slope: f64;
        let right_slope: f64;

        if self.config.using_numerov {
            let left_last_index = self.left_wavefunction.len() - 1;
            let psi_left_one_next = self.next(
                &Side::Left,
                left_last_index,
                self.left_wavefunction[left_last_index],
                self.left_wavefunction[left_last_index - 1],
            );
            let psi_left_two_next = self.next(
                &Side::Left,
                left_last_index + 1,
                psi_left_one_next,
                self.left_wavefunction[left_last_index],
            );

            let psi_left_one_prev = self.left_wavefunction[left_last_index - 1];
            let psi_left_two_prev = self.left_wavefunction[left_last_index - 2];

            left_slope = (-psi_left_two_next + 8.0 * psi_left_one_next - 8.0 * psi_left_one_prev
                + psi_left_two_prev)
                / (12.0 * self.config.step_size);

            let right_last_index = self.right_wavefunction.len() - 1;
            let psi_right_one_next = self.next(
                &Side::Right,
                right_last_index,
                self.right_wavefunction[right_last_index],
                self.right_wavefunction[right_last_index - 1],
            );
            let psi_right_two_next = self.next(
                &Side::Right,
                right_last_index + 1,
                psi_right_one_next,
                self.right_wavefunction[right_last_index],
            );

            let psi_right_one_prev = self.right_wavefunction[right_last_index - 1];
            let psi_right_two_prev = self.right_wavefunction[right_last_index - 2];

            right_slope = (-psi_right_two_prev + 8.0 * psi_right_one_prev
                - 8.0 * psi_right_one_next
                + psi_right_two_next)
                / (12.0 * self.config.step_size);
        } else {
            left_slope = self.left_wavefunction.last().unwrap()
                - self.left_wavefunction[self.left_wavefunction.len() - 2];

            right_slope = -(self.right_wavefunction.last().unwrap()
                - self.right_wavefunction[self.right_wavefunction.len() - 2]);
        }
        (left_slope, right_slope)
    }

    /// Normalizes the wavefunction so that ∫|ψ|²dx = 1.
    fn normalize(&mut self) {
        let (_, mut f): (Vec<_>, Vec<_>) = self.wavefunction_points().iter().cloned().unzip();
        f = f.iter().map(|val| val * val).collect();
        let integral = trapezoidal(&f, &self.config.step_size);
        self.left_wavefunction
            .iter_mut()
            .for_each(|val| *val = *val * (1.0 / integral).sqrt());

        self.right_wavefunction
            .iter_mut()
            .for_each(|val| *val = *val * (1.0 / integral).sqrt());
    }
}

impl Solver for MatchingSolver {
    type CONFIG = MatchingConfig;

    fn new(config: &Self::CONFIG) -> Self {
        let steps = ((config.x_max - config.x_min) / config.step_size).round() as usize + 1;
        let match_idx = ((config.x_match - config.x_min) / config.step_size).round() as usize;

        MatchingSolver {
            config: config.clone(),
            steps,
            energy: config.initial_energy,
            energy_step_size: config.initial_energy_step_size,
            is_left_slope_larger: None,
            left_wavefunction: Vec::<f64>::with_capacity(match_idx + 1),
            right_wavefunction: Vec::<f64>::with_capacity(steps - match_idx),
        }
    }

    /// Popuplates the wavefunction vector with a solution to the Schrodinger equation
    /// and also determines the corresponding energy. The process requires iterating
    /// over many candidate energies and stopping when the energy step size becomes
    /// sufficiently small.
    fn solve(&mut self) {
        loop {
            let result = self.compute_wavefunction();
            let (_, f): (Vec<_>, Vec<_>) = self.wavefunction_points().iter().cloned().unzip();
            //let e = energy_of(&f, self.steps, self.step_size, self.potential, self.x_min);
            //dbg!(self.energy, e);

            if self.energy_step_size.abs() <= self.config.energy_step_size_cutoff {
                self.normalize();
                break;
            }

            let (left_slope, right_slope) = self.slopes();

            if result.is_ok() && self.is_left_slope_larger == Some(left_slope < right_slope) {
                self.energy_step_size = -self.energy_step_size / 2.0;
            }

            self.is_left_slope_larger = Some(left_slope >= right_slope);

            self.energy += self.energy_step_size;
        }
    }

    fn energy(&self) -> f64 {
        self.energy
    }

    fn reset(&mut self) {
        *self = Self::new(&self.config)
    }

    fn wavefunction_points(&self) -> Vec<(f64, f64)> {
        let mut pairs: Vec<(f64, f64)> = Vec::with_capacity(self.steps);

        for (i, psi_val) in self.left_wavefunction.iter().enumerate() {
            pairs.push((self.x_from_index(i, &Side::Left), *psi_val));
        }

        for (i, psi_val) in self.right_wavefunction.iter().enumerate().rev().skip(1) {
            pairs.push((self.x_from_index(i, &Side::Right), *psi_val));
        }

        pairs
    }
}

/// Indicates the side of a component wavefunction in the matching method.
enum Side {
    Left,
    Right,
}

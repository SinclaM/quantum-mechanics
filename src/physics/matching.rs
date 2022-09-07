//! Matching method for solving the time-independent Schrodinger equation
//! in one dimension.

use crate::utils::integration::trapezoidal;
use std::io::Write;

/// A solver that looks for solutions using the matching method.

pub struct MatchingSolver {
    pub steps: usize,
    pub step_size: f64,
    pub energy: f64,
    pub energy_step_size: f64,
    pub potential: fn(f64) -> f64,
    pub energy_step_size_cutoff: f64,
    is_left_slope_larger: bool,
    pub left_wavefunction: Vec<f64>,
    pub right_wavefunction: Vec<f64>,
    match_idx: usize,
    pub x_min: f64,
    pub x_max: f64,
    pub using_numerov: bool,
}

impl MatchingSolver {
    /// Creates a new matching solver.
    pub fn new(
        step_size: f64,
        energy: f64,
        energy_step_size: f64,
        potential: fn(f64) -> f64,
        energy_step_size_cutoff: f64,
        x_min: f64,
        x_max: f64,
        match_idx: usize,
        using_numerov: bool,
    ) -> MatchingSolver {
        let steps = ((x_max - x_min) / step_size).round() as usize + 1;

        MatchingSolver {
            steps,
            step_size,
            energy,
            energy_step_size,
            potential,
            energy_step_size_cutoff,
            is_left_slope_larger: true,
            left_wavefunction: Vec::<f64>::with_capacity(match_idx + 1),
            right_wavefunction: Vec::<f64>::with_capacity(steps - match_idx),
            match_idx,
            x_min,
            x_max,
            using_numerov,
        }
    }

    fn x_from_index(&self, i: usize, side: &Side) -> f64 {
        match side {
            Side::Left  => self.x_min + (i as f64) * self.step_size,
            Side::Right => self.x_max - (i as f64) * self.step_size
        }
    }

    /// Computes a term needed for the Numerov method.
    fn k_sqr(&self, x: f64) -> f64 {
        2.0 * (self.energy - (self.potential)(x))
    }

    /// Applies the finite difference approximation to find the value of wavefunction
    /// one position toward the matching point from either the left or right.
    fn step(&mut self, side: &Side) {
        let next: f64;
        match side {
            Side::Left => {
                let last_index = self.left_wavefunction.len() - 1;
                if self.using_numerov {
                    next = (2.0
                        * (1.0
                            - (5.0 / 12.0)
                                * self.step_size.powf(2.0)
                                * self.k_sqr(self.x_from_index(last_index, &side)))
                        * self.left_wavefunction[last_index]
                        - (1.0
                            + (1.0 / 12.0)
                                * self.step_size.powf(2.0)
                                * self.k_sqr(self.x_from_index(last_index - 1, &side)))
                            * self.left_wavefunction[last_index - 1])
                        / (1.0
                            + (1.0 / 12.0)
                                * self.step_size.powf(2.0)
                                * self.k_sqr(self.x_from_index(last_index + 1, &side)));
                } else {
                    next = 2.0
                        * (self.step_size
                            * self.step_size
                            * ((self.potential)(self.x_from_index(last_index, &side)) - self.energy)
                            + 1.0)
                        * self.left_wavefunction[last_index]
                        - self.left_wavefunction[last_index - 1];
                }
                self.left_wavefunction.push(next);
            },
            Side::Right => {
                let last_index = self.right_wavefunction.len() - 1;
                if self.using_numerov {
                    next = (2.0
                        * (1.0
                            - (5.0 / 12.0)
                                * self.step_size.powf(2.0)
                                * self.k_sqr(self.x_from_index(last_index, &side)))
                        * self.right_wavefunction[last_index]
                        - (1.0
                            + (1.0 / 12.0)
                                * self.step_size.powf(2.0)
                                * self.k_sqr(self.x_from_index(last_index - 1, &side)))
                            * self.right_wavefunction[last_index - 1])
                        / (1.0
                            + (1.0 / 12.0)
                                * self.step_size.powf(2.0)
                                * self.k_sqr(self.x_from_index(last_index + 1, &side)));
                } else {
                    next = 2.0
                        * (self.step_size
                            * self.step_size
                            * ((self.potential)(self.x_from_index(last_index, &side)) - self.energy)
                            + 1.0)
                        * self.right_wavefunction[last_index]
                        - self.right_wavefunction[last_index - 1];
                }
                self.right_wavefunction.push(next);
            }
        }

    }

    /// Approximates the wavefunction for the current energy. Stops when it has
    /// computed the requested number of steps, or if the wavefunction begins
    /// diverging.
    fn compute_wavefunction(&mut self) {
        self.reset_wavefunction();
        for _ in 1..=(self.match_idx - 1) {
            self.step(&Side::Left);
        }

        for _ in 1..=(self.steps - self.match_idx - 2) {
            self.step(&Side::Right);
        }

        // Ensure continuity at match_idx
        let scale_factor =
            self.left_wavefunction.last().unwrap() / self.right_wavefunction.last().unwrap();
        self.right_wavefunction
            .iter_mut()
            .for_each(|val| *val *= scale_factor);
    }

    /// Resets the wavefunction vectors.
    fn reset_wavefunction(&mut self) {
        self.left_wavefunction.clear();
        self.right_wavefunction.clear();

        self.left_wavefunction.push(0.0);
        self.left_wavefunction.push(0.001 * self.step_size);

        self.right_wavefunction.push(0.0);
        self.right_wavefunction.push(0.001 * self.step_size);
    }

    /// Popuplates the wavefunction vector with a solution to the Schrodinger equation
    /// and also determines the corresponding energy. The process requires iterating
    /// over many candidate energies and stopping when the energy step size becomes
    /// sufficiently small.
    pub fn solve(&mut self) {
        loop {
            self.compute_wavefunction();
            if self.energy_step_size.abs() <= self.energy_step_size_cutoff {
                self.normalize();
                break;
            }

            let left_rise = self.left_wavefunction.last().unwrap()
                - self.left_wavefunction[self.left_wavefunction.len() - 2];

            let right_rise = -(self.right_wavefunction.last().unwrap()
                - self.right_wavefunction[self.right_wavefunction.len() - 2]);

            if self.is_left_slope_larger != (left_rise >= right_rise) {
                self.energy_step_size = -self.energy_step_size / 2.0;
            }

            self.is_left_slope_larger = left_rise >= right_rise;

            self.energy += self.energy_step_size;
        }
    }

    /// Prints energy and wavefunction data to a text file. Useful for later analysis
    /// with tools like gnuplot. The first line is the a `#` (gnuplot comment) followed
    /// by the energy value (e.g "# 1.24345678"). The rest of the lines in the file
    /// are the x-value and then the wavefunction value (separated by a space).
    pub fn dump_to_file(&mut self, data_file: &mut std::fs::File) -> Result<(), std::io::Error> {
        writeln!(data_file, "# {}", self.energy)?;

        for (x, psi) in self.wavefunction_points().iter() {
            writeln!(data_file, "{} {}", x, psi)?;
        }

        Ok(())
    }

    /// Returns a vector of (x, ψ) points for the wavefunction.
    pub fn wavefunction_points(&self) -> Vec<(f64, f64)> {
        let mut pairs: Vec<(f64, f64)> = Vec::with_capacity(self.steps);
        for (i, psi_val) in self.left_wavefunction.iter().enumerate() {
            pairs.push((self.x_from_index(i, &Side::Left), *psi_val));
        }

        for (i, psi_val) in self.right_wavefunction.iter().enumerate().rev().skip(1) {
            pairs.push((self.x_from_index(i, &Side::Right), *psi_val));
        }

        pairs
    }

    /// Normalizes the wavefunction so that ∫|ψ|²dx = 1.
    fn normalize(&mut self) {
        let (_, mut f): (Vec<_>, Vec<_>) = self.wavefunction_points().iter().cloned().unzip();
        f = f.iter().map(|val| val * val).collect();
        let integral = trapezoidal(&f, &self.step_size);
        self.left_wavefunction
            .iter_mut()
            .for_each(|val| *val = *val * (1.0 / integral).sqrt());

        self.right_wavefunction
            .iter_mut()
            .for_each(|val| *val = *val * (1.0 / integral).sqrt());
    }
}

enum Side {
    Left,
    Right
}

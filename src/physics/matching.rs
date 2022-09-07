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
    using_numerov: bool,
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

    /// Computes the x value associated with an index into the left wavefunction.
    fn x_from_left_idx(&self, i: usize) -> f64 {
        self.x_min + (i as f64) * self.step_size
    }

    /// Computes the x value associated with an index into the right wavefunction.
    fn x_from_right_idx(&self, i: usize) -> f64 {
        self.x_max - (i as f64) * self.step_size
    }

    /// Applies the finite difference approximation to find the value of wavefunction
    /// one position toward the matching point from the left.
    fn step_left(&mut self) {
        let last_index = self.left_wavefunction.len() - 1;
        let mut next: f64;
        if self.using_numerov {
            next =
                2.0 * self.numerov_factor(
                    self.x_from_left_idx(last_index),
                    self.left_wavefunction[last_index],
                ) - self.numerov_factor(
                    self.x_from_left_idx(last_index - 1),
                    self.left_wavefunction[last_index - 1],
                ) - 2.0
                    * (self.energy - (self.potential)(self.x_from_left_idx(last_index)))
                    * (self.step_size * self.step_size)
                    * self.numerov_factor(
                        self.x_from_left_idx(last_index),
                        self.left_wavefunction[last_index],
                    );
            // Need to recover ψ from Y.
            next /= 1.0
                - (1.0 / 6.0)
                    * (self.step_size * self.step_size)
                    * ((self.potential)(self.x_from_left_idx(last_index + 1)) - self.energy);
        } else {
            next = 2.0 * self.left_wavefunction[last_index]
                - self.left_wavefunction[last_index - 1]
                - 2.0
                    * (self.energy - (self.potential)(self.x_from_left_idx(last_index)))
                    * (self.step_size * self.step_size)
                    * self.left_wavefunction[last_index];
        }

        self.left_wavefunction.push(next);
    }

    /// Applies the finite difference approximation to find the value of wavefunction
    /// one position toward the matching point from the right.
    fn step_right(&mut self) {
        let last_index = self.right_wavefunction.len() - 1;
        let mut next: f64;
        if self.using_numerov {
            next = 2.0
                * (self.step_size
                    * self.step_size
                    * ((self.potential)(self.x_from_right_idx(last_index)) - self.energy)
                    + 1.0)
                * self.numerov_factor(
                    self.x_from_right_idx(last_index),
                    self.right_wavefunction[last_index],
                )
                - self.numerov_factor(
                    self.x_from_right_idx(last_index - 1),
                    self.right_wavefunction[last_index - 1],
                );
            // Need to recover ψ from Y.
            next /= 1.0
                - (1.0 / 6.0)
                    * (self.step_size * self.step_size)
                    * ((self.potential)(self.x_from_right_idx(last_index + 1)) - self.energy);
        } else {
            next = 2.0
                * (self.step_size
                    * self.step_size
                    * ((self.potential)(self.x_from_right_idx(last_index)) - self.energy)
                    + 1.0)
                * self.right_wavefunction[last_index]
                - self.right_wavefunction[last_index - 1];
        }
        self.right_wavefunction.push(next);
    }

    /// Computes the numerov factor, Y.
    fn numerov_factor(&self, x: f64, psi: f64) -> f64 {
        (1.0 - (1.0 / 6.0)
            * (self.step_size * self.step_size)
            * ((self.potential)(x) - self.energy))
            * psi
    }

    /// Approximates the wavefunction for the current energy. Stops when it has
    /// computed the requested number of steps, or if the wavefunction begins
    /// diverging.
    fn compute_wavefunction(&mut self) {
        self.reset_wavefunction();
        for _ in 1..=(self.match_idx - 1) {
            self.step_left();
        }

        for _ in 1..=(self.steps - self.match_idx - 2) {
            self.step_right();
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
            pairs.push((self.x_from_left_idx(i), *psi_val));
        }

        for (i, psi_val) in self.right_wavefunction.iter().enumerate().rev().skip(1) {
            pairs.push((self.x_from_right_idx(i), *psi_val));
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

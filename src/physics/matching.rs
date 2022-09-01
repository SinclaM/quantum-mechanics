//! Matching method for solving the time-independent Schrodinger equation
//! in one dimension.

use crate::utils::{InwardVec, InwardVecError};
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
    pub wavefunction: InwardVec<f64>,
    match_idx: usize,
    pub x_min: f64,
    pub x_max: f64,
}

impl MatchingSolver {
    /// Returns a new solver for even parity wavefunctions.
    pub fn new(
        step_size: f64,
        energy: f64,
        energy_step_size: f64,
        potential: fn(f64) -> f64,
        energy_step_size_cutoff: f64,
        x_min: f64,
        x_max: f64,
        match_idx: usize,
    ) -> Result<MatchingSolver, InwardVecError> {
        let steps = ((x_max - x_min) / step_size) as usize;

        Ok(MatchingSolver {
            steps,
            step_size,
            energy,
            energy_step_size,
            potential,
            energy_step_size_cutoff,
            is_left_slope_larger: true,
            wavefunction: InwardVec::<f64>::new(steps, match_idx)?,
            match_idx,
            x_min,
            x_max,
        })
    }

    /// Applies the finite difference approximation to find value of wavefunction
    /// one position toward the matching point from the left.
    fn step_from_left(&mut self) -> Result<(), InwardVecError> {
        let last_left_idx = self.wavefunction.left - 1;
        self.wavefunction.push_from_left(
            2.0 * self.wavefunction.data[last_left_idx]
                - self.wavefunction.data[last_left_idx - 1]
                - 2.0
                    * (self.energy
                        - (self.potential)(self.x_min + (last_left_idx as f64) * self.step_size))
                    * (self.step_size * self.step_size)
                    * self.wavefunction.data[last_left_idx],
        )?;
        Ok(())
    }

    /// Applies the finite difference approximation to find value of wavefunction
    /// one position toward the matching point from the right.
    fn step_from_right(&mut self) -> Result<(), InwardVecError> {
        let last_right_index = self.wavefunction.right + 1;
        self.wavefunction.push_from_right(
            2.0 * (self.step_size
                * self.step_size
                * ((self.potential)(self.x_min + (last_right_index as f64) * self.step_size) - self.energy)
                + 1.0)
                * self.wavefunction.data[last_right_index]
                - self.wavefunction.data[last_right_index + 1],
        )?;
        Ok(())
    }

    /// Approximates the wavefunction for the current energy. Stops when it has
    /// computed the requested number of steps, or if the wavefunction begins
    /// diverging.
    fn compute_wavefunction(&mut self) {
        self.reset_wavefunction();

        while self.step_from_left().is_ok() {
            continue;
        }

        while self.step_from_right().is_ok() {
            continue;
        }

        // Ensure continuity at match_idx
        let last_right_index = self.wavefunction.right + 1;
        let right_at_match = 2.0
            * (self.step_size
                * self.step_size
                * ((self.potential)(self.x_min + (last_right_index as f64) * self.step_size)
                    - self.energy)
                + 1.0)
            * self.wavefunction.data[last_right_index]
            - self.wavefunction.data[last_right_index + 1];

        let scale_factor = self.wavefunction.data[self.wavefunction.left - 1] / right_at_match;
        self.wavefunction
            .data
            .iter_mut()
            .skip(self.match_idx + 1)
            .for_each(|val| *val *= scale_factor);
    }

    /// Resets the wavefunction vector. The values at the first two positions are
    /// set to 1.0, which is appopriate for non-normalized even parity wavefunctions.
    fn reset_wavefunction(&mut self) {
        self.wavefunction
            .data
            .iter_mut()
            .for_each(|val| *val = Default::default());
        self.wavefunction.left = 0;
        self.wavefunction.right = self.wavefunction.data.len() - 1;

        self.wavefunction.push_from_left(0.0).unwrap();
        self.wavefunction
            .push_from_left(0.0001 * self.step_size)
            .unwrap();

        self.wavefunction.push_from_right(0.0).unwrap();
        self.wavefunction
            .push_from_right(0.0001 * self.step_size)
            .unwrap();
    }

    /// Popuplates the wavefunction vector with a solution to the Schrodinger equation
    /// and also determines the corresponding energy. The process requires iterating
    /// over many candidate energies and stopping when the energy step size becomes
    /// sufficiently small.
    pub fn solve(&mut self) {
        loop {
            self.compute_wavefunction();
            if self.energy_step_size.abs() <= self.energy_step_size_cutoff {
                break;
            }

            let left_slope = self.wavefunction.data[self.wavefunction.left - 1]
                - self.wavefunction.data[self.wavefunction.left - 2];

            let last_right_index = self.wavefunction.right + 1;
            let right_at_match = 2.0
                * (self.step_size
                    * self.step_size
                    * ((self.potential)(self.x_min + (last_right_index as f64) * self.step_size)
                        - self.energy)
                    + 1.0)
                * self.wavefunction.data[last_right_index]
                - self.wavefunction.data[last_right_index + 1];

            let right_slope =
                -(right_at_match - self.wavefunction.data[self.wavefunction.right + 1]);

            if self.is_left_slope_larger != (left_slope >= right_slope) {
                self.energy_step_size = -self.energy_step_size / 2.0;
            }

            self.is_left_slope_larger = left_slope >= right_slope;

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
        for (i, psi_val) in self.wavefunction.data.iter().enumerate() {
            pairs.push((self.x_min + (i as f64) * self.step_size, *psi_val));
        }

        pairs
    }
}

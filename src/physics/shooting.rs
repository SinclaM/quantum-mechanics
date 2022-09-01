//! Shooting method for solving the time-independent Schrodinger equation
//! in one dimension. The shooting method is only applicable for even
//! potentials (symmetric about x = 0).

use std::io::Write;

/// A solver that looks for solutions of the desired parity
/// using the shooting method.

pub struct ShootingSolver {
    pub steps: usize,
    pub step_size: f64,
    pub energy: f64,
    pub energy_step_size: f64,
    pub wavefunction_cutoff: f64,
    pub wavefunction: Vec<f64>,
    last_diverge: f64,
    pub potential: fn(f64) -> f64,
    pub energy_step_size_cutoff: f64,
    pub parity: Parity,
}

impl ShootingSolver {
    /// Returns a new solver for even parity wavefunctions.
    pub fn new(
        steps: usize,
        step_size: f64,
        energy: f64,
        energy_step_size: f64,
        wavefunction_cutoff: f64,
        potential: fn(f64) -> f64,
        energy_step_size_cutoff: f64,
        parity: Parity,
    ) -> ShootingSolver {
        let wavefunction: Vec<f64> = Vec::with_capacity(steps);

        ShootingSolver {
            steps,
            step_size,
            energy,
            energy_step_size,
            wavefunction_cutoff,
            wavefunction,
            last_diverge: 0.0,
            potential,
            energy_step_size_cutoff,
            parity,
        }
    }

    /// Returns a new solver for even parity wavefunctions, using default options
    /// for some fields.
    pub fn default(
        steps: usize,
        step_size: f64,
        energy: f64,
        potential: fn(f64) -> f64,
        parity: Parity,
    ) -> ShootingSolver {
        let wavefunction: Vec<f64> = Vec::with_capacity(10000);

        ShootingSolver {
            steps,
            step_size,
            energy,
            energy_step_size: 0.1,
            wavefunction_cutoff: 100.0,
            wavefunction,
            last_diverge: 0.0,
            potential,
            energy_step_size_cutoff: 0.000001,
            parity,
        }
    }

    /// Applies the finite difference approximation to find value of wavefunction
    /// one position forward, using the two most recent values.
    fn step(&mut self) {
        let i = self.wavefunction.len() - 1;
        self.wavefunction.push(
            2.0 * self.wavefunction[i]
                - self.wavefunction[i - 1]
                - 2.0
                    * (self.energy - (self.potential)((i as f64) * self.step_size))
                    * (self.step_size * self.step_size)
                    * self.wavefunction[i],
        );
    }

    /// Determines if the wavefunction is diverging to infinity (positive or negative).
    fn is_diverging(&self) -> bool {
        self.wavefunction.last().unwrap().abs() > self.wavefunction_cutoff
    }

    /// Approximates the wavefunction for the current energy. Stops when it has
    /// computed the requested number of steps, or if the wavefunction begins
    /// diverging.
    fn compute_wavefunction(&mut self) {
        self.reset_wavefunction();
        for _ in 2..=self.steps {
            if self.is_diverging() {
                break;
            }
            self.step();
        }
    }

    /// Resets the wavefunction vector. The values at the first two positions are
    /// set to 1.0, which is appopriate for non-normalized even parity wavefunctions.
    fn reset_wavefunction(&mut self) {
        self.wavefunction.clear();
        match self.parity {
            Parity::Even => {
                self.wavefunction.push(1.0);
                self.wavefunction.push(1.0);
            }
            Parity::Odd => {
                self.wavefunction.push(0.0);
                self.wavefunction.push(self.step_size);
            }
        }
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

            if self.wavefunction.last().unwrap() * self.last_diverge < 0.0 {
                self.energy_step_size = -self.energy_step_size / 2.0;
            }

            self.energy += self.energy_step_size;
            self.last_diverge = if self.wavefunction.last().unwrap() >= &0.0 {
                1.0
            } else {
                -1.0
            };
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
        let mut x_vals: Vec<f64> = Vec::new();
        let mut psi_vals: Vec<f64> = Vec::new();
        let mut pairs: Vec<(f64, f64)> = Vec::new();

        x_vals.push(-(self.wavefunction.len() as f64) * self.step_size);

        for val in self.wavefunction.iter().skip(1).rev() {
            x_vals.push(x_vals.last().unwrap() + self.step_size);
            psi_vals.push(match self.parity {
                Parity::Even => *val,
                Parity::Odd => -val,
            });
            pairs.push((*x_vals.last().unwrap(), *psi_vals.last().unwrap()));
        }
        for val in &self.wavefunction {
            x_vals.push(x_vals.last().unwrap() + self.step_size);
            psi_vals.push(*val);
            pairs.push((*x_vals.last().unwrap(), *psi_vals.last().unwrap()));
        }
        pairs
    }
}

/// The parity of solutions that a solver will look for when solving
/// the Schrodinger equation.
pub enum Parity {
    Even,
    Odd,
}

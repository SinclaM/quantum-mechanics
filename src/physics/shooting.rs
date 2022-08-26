//! Shooting method for solving the time-independent Schrodinger equation
//! in one dimension.

/// A solver that looks for even parity solutions using the shooting
/// method.
pub struct ShootingSolverEven {
    pub steps: usize,
    pub step_size: f64,
    pub energy: f64,
    pub energy_step_size: f64,
    pub wavefunction_cutoff: f64,
    pub wavefunction: Vec<f64>,
    last_diverge: f64,
    pub potential: fn(f64) -> f64,
    pub energy_step_size_cutoff: f64,
}

impl ShootingSolverEven {
    /// Returns a new solver for even parity wavefunctions.
    pub fn new(
        steps: usize,
        step_size: f64,
        energy: f64,
        energy_step_size: f64,
        wavefunction_cutoff: f64,
        potential: fn(f64) -> f64,
        energy_step_size_cutoff: f64,
    ) -> ShootingSolverEven {
        let wavefunction: Vec<f64> = Vec::with_capacity(steps);

        ShootingSolverEven {
            steps,
            step_size,
            energy,
            energy_step_size,
            wavefunction_cutoff,
            wavefunction,
            last_diverge: 0.0,
            potential,
            energy_step_size_cutoff,
        }
    }

    /// Returns a new solver for even parity wavefunctions, using default options
    /// for some fields.
    pub fn default(
        steps: usize,
        step_size: f64,
        energy: f64,
        potential: fn(f64) -> f64,
    ) -> ShootingSolverEven {
        let wavefunction: Vec<f64> = Vec::with_capacity(10000);

        ShootingSolverEven {
            steps,
            step_size,
            energy,
            energy_step_size: 10.0,
            wavefunction_cutoff: 300.0,
            wavefunction,
            last_diverge: 0.0,
            potential,
            energy_step_size_cutoff: 0.000001,
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
        for _ in 2..self.steps {
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
        self.wavefunction.push(1.0);
        self.wavefunction.push(1.0);
    }

    /// Popuplates the wavefunction vector with a solution to the Schrodinger equation
    /// and also determines the corresponding energy. The process requires iterating
    /// over many candidate energies and stopping when the energy step size becomes
    /// sufficiently small.
    pub fn solve(&mut self) {
        loop {
            self.reset_wavefunction();
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
}

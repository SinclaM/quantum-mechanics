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

    fn is_diverging(&self) -> bool {
        self.wavefunction.last().unwrap().abs() > self.wavefunction_cutoff
    }

    fn compute_wavefunction(&mut self) {
        for _ in 2..self.steps {
            if self.is_diverging() {
                break;
            }
            self.step();
        }
    }

    fn reset_wavefunction(&mut self) {
        self.wavefunction.clear();
        self.wavefunction.push(1.0);
        self.wavefunction.push(1.0);
    }

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

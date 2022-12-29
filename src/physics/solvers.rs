use std::io::Write;

pub mod shooting;
pub mod matching;
pub mod variational;

pub trait Solver {
    type CONFIG;

    fn new(config: &Self::CONFIG) -> Self;

    fn solve(&mut self);

    fn energy(&self) -> f64;

    /// Resets the solver to it's initial configuration
    fn reset(&mut self);

    /// Returns a vector of (x, ψ) points for the wavefunction.
    fn wavefunction_points(&self) -> Vec<(f64, f64)>;

    /// Prints energy and wavefunction data to a text file. Useful for later analysis
    /// with tools like gnuplot. The first line is the a '#' (gnuplot comment) followed
    /// by the energy value (e.g "# 1.24345678"). The rest of the lines in the file
    /// are the x-value and then the wavefunction value (separated by a space).
    /// # Example
    /// ```txt
    /// --output file--
    ///     # 22.0927734375
    ///     -1.3 0
    ///     -1.25 0.00000000928011292465171
    ///     -1.2 -0.00000009554014322853878
    ///     -1.15 0.00000097431987580849
    ///     -1.1 -0.000009935227934806336
    ///       ︙             ︙
    /// ```
    fn dump_to_file(&mut self, data_file: &mut std::fs::File) -> Result<(), std::io::Error> {
        writeln!(data_file, "# {}", self.energy())?;

        for (x, psi) in self.wavefunction_points().iter() {
            writeln!(data_file, "{} {}", x, psi)?;
        }

        Ok(())
    }
}

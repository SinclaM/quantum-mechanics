pub mod integration;
pub mod root_finding;


pub fn gen_range(range: std::ops::RangeInclusive<f64>, step: f64) -> Vec<f64> {
    let mut x = *range.start();
    let mut ret = Vec::<f64>::new();
    while x <= *range.end() {
        ret.push(x);
        x += step;
    }
    ret
}

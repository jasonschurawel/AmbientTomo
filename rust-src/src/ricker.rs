pub fn ricker_wavelet(t: f64, f_dominant: f64) -> f64 {
    let t0 = 1.2 / f_dominant;
    let tau = std::f64::consts::PI * f_dominant * (t - t0);
    (1.0 - 2.0 * tau * tau) * (-tau * tau).exp()
}
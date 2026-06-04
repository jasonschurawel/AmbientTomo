pub fn erzeuge_cpml_koeffizienten(
    n: usize,
    pml_punkte: usize,
    d: f64,
    v_max: f64,
    dt: f64,
    f_dominant: f64,
) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let r_theor = 1e-6f64;
    let l = (pml_punkte as f64) * d;
    let sigma_max = -(3.0 * v_max * r_theor.ln()) / (2.0 * l);
    let alpha_max = std::f64::consts::PI * f_dominant;

    let mut sigma_int = vec![0.0; n.saturating_sub(2)];
    let mut alpha_int = vec![0.0; n.saturating_sub(2)];

    for idx in 0..n.saturating_sub(2) {
        let g_idx = idx + 1;
        let dist = if g_idx < pml_punkte {
            ((pml_punkte - g_idx) as f64) * d
        } else if g_idx >= n - pml_punkte {
            ((g_idx - (n - 1 - pml_punkte)) as f64) * d
        } else {
            0.0
        };
        sigma_int[idx] = sigma_max * (dist / l).powi(2);
        alpha_int[idx] = alpha_max * (1.0 - dist / l);
    }

    let mut sigma_half = vec![0.0; n.saturating_sub(1)];
    let mut alpha_half = vec![0.0; n.saturating_sub(1)];

    for idx in 0..n.saturating_sub(1) {
        let h_idx = (idx as f64) + 0.5;
        let dist = if h_idx < pml_punkte as f64 {
            ((pml_punkte as f64) - h_idx) * d
        } else if h_idx > (n - 1 - pml_punkte) as f64 {
            (h_idx - (n - 1 - pml_punkte) as f64) * d
        } else {
            0.0
        };
        sigma_half[idx] = sigma_max * (dist / l).powi(2);
        alpha_half[idx] = alpha_max * (1.0 - dist / l);
    }

    let b_int: Vec<f64> = sigma_int.iter().zip(&alpha_int)
        .map(|(&s, &a)| (-(s + a) * dt).exp()).collect();
    let mut a_int = vec![0.0; n.saturating_sub(2)];
    for i in 0..a_int.len() {
        let sum = sigma_int[i] + alpha_int[i];
        if sum > 0.0 {
            a_int[i] = (b_int[i] - 1.0) * sigma_int[i] / sum;
        }
    }

    let b_half: Vec<f64> = sigma_half.iter().zip(&alpha_half)
        .map(|(&s, &a)| (-(s + a) * dt).exp()).collect();
    let mut a_half = vec![0.0; n.saturating_sub(1)];
    for i in 0..a_half.len() {
        let sum = sigma_half[i] + alpha_half[i];
        if sum > 0.0 {
            a_half[i] = (b_half[i] - 1.0) * sigma_half[i] / sum;
        }
    }

    (a_int, b_int, a_half, b_half)
}
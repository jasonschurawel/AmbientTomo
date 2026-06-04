#[derive(Clone)]
pub struct Params {
    pub t_sim_sec: f64,
    pub breite_m: f64,
    pub tiefe_m: f64,
    pub grenzschicht_m: f64,
    pub pml_dicke_m: f64,
    pub quelle_x_m: f64,
    pub quelle_tiefe_m: f64,
    pub v_oben: f64,
    pub v_unten: f64,
    pub f_dominant: f64,
    pub f_max_effektiv: f64,
    pub punkte_pro_welle: f64,
    pub dampfung_global: f64,
}
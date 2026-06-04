use std::env;
use std::process;

mod cpml;
mod params;
mod ricker;
mod simulation;
mod render;

use params::Params;

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: {} <output.mp4>", args[0]);
        process::exit(1);
    }

    let output_video = &args[1];

    let params = Params {
        t_sim_sec: 2.5,
        breite_m: 5000.0,
        tiefe_m: 1000.0,
        grenzschicht_m: 500.0,
        pml_dicke_m: 250.0,
        quelle_x_m: 400.0,
        quelle_tiefe_m: 100.0,
        v_oben: 1500.0,
        v_unten: 3000.0,
        f_dominant: 15.0,
        f_max_effektiv: 40.0,
        punkte_pro_welle: 10.0,
        dampfung_global: 0.4,
    };

    println!("🚀 Starting Rust seismic simulation...");

    if let Err(e) = simulation::wellen_simulieren_2d_cpml(&params, "temp.bin") {
        eprintln!("Simulation error: {}", e);
        process::exit(1);
    }

    println!("🎥 Rendering video...");
    if let Err(e) = render::create_video("temp.bin", output_video, &params) {
        eprintln!("Rendering error: {}", e);
        process::exit(1);
    }

    let _ = std::fs::remove_file("temp.bin");
    println!("✅ Done! Video saved as: {}", output_video);
}
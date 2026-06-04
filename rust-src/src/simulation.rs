use crate::cpml::erzeuge_cpml_koeffizienten;
use crate::ricker::ricker_wavelet;
use crate::Params;
use ndarray::{Array2, s, Axis};
use rayon::prelude::*;
use std::fs::File;
use std::io::{self, Write};

pub fn wellen_simulieren_2d_cpml(params: &Params, bin_path: &str) -> io::Result<()> {
    let lambda_min = params.v_oben / params.f_max_effektiv;
    let dx = lambda_min / params.punkte_pro_welle;
    let dy = dx;

    let pml = (params.pml_dicke_m / dx) as usize;
    let nx = (params.breite_m / dx) as usize + 2 * pml;
    let ny = (params.tiefe_m / dy) as usize + 2 * pml;

    let mut v = Array2::from_elem((nx, ny), params.v_oben);
    let grenz = (params.grenzschicht_m / dy) as usize + pml;
    v.slice_mut(s![.., grenz..]).fill(params.v_unten);

    let dt = 0.5 / (params.v_unten * (1.0/dx.powi(2) + 1.0/dy.powi(2)).sqrt());
    let nt = (params.t_sim_sec / dt) as usize;

    let qx = (params.quelle_x_m / dx) as usize + pml;
    let qy = (params.quelle_tiefe_m / dy) as usize + pml;

    println!("Grid: {}x{} | dt={:.6}s | PML={}", nx, ny, dt, pml);

    let (a_x_int, b_x_int, a_x_half, b_x_half) = 
        erzeuge_cpml_koeffizienten(nx, pml, dx, params.v_unten, dt, params.f_dominant);
    let (a_y_int, b_y_int, a_y_half, b_y_half) = 
        erzeuge_cpml_koeffizienten(ny, pml, dy, params.v_unten, dt, params.f_dominant);

    let mut u = Array2::<f64>::zeros((nx, ny));
    let mut u_past = Array2::<f64>::zeros((nx, ny));
    let mut u_next = Array2::<f64>::zeros((nx, ny));

    let mut psi_x = Array2::<f64>::zeros((nx-1, ny-2));
    let mut phi_x = Array2::<f64>::zeros((nx-2, ny-2));
    let mut psi_y = Array2::<f64>::zeros((nx-2, ny-1));
    let mut phi_y = Array2::<f64>::zeros((nx-2, ny-2));

    let mut file = File::create(bin_path)?;

    file.write_all(&dx.to_le_bytes())?;
    file.write_all(&dy.to_le_bytes())?;
    file.write_all(&dt.to_le_bytes())?;
    file.write_all(&(pml as u32).to_le_bytes())?;
    file.write_all(&(nx as u32).to_le_bytes())?;
    file.write_all(&(ny as u32).to_le_bytes())?;

    let save_step = (nt / 900).max(1);

    for i in 0..nt.saturating_sub(2) {
        if i % (nt / 20).max(1) == 0 {
            print!("\rSim: [{:3.0}%]", (i as f64 / nt as f64) * 100.0);
            let _ = std::io::stdout().flush();
        }

        // === CPML X Direction (Berechnung der gestreckten Ableitung) ===
        let d2u_dx2_stretched = {
            let du_dx = (&u.slice(s![1.., 1..-1]) - &u.slice(s![..-1, 1..-1])) / dx;
            let b2d = Array2::from_shape_vec((nx-1, 1), b_x_half.clone()).unwrap();
            let a2d = Array2::from_shape_vec((nx-1, 1), a_x_half.clone()).unwrap();
            psi_x.assign(&(&b2d * &psi_x + &a2d * &du_dx));

            let hx = &du_dx + &psi_x;
            let dh_dx = (&hx.slice(s![1.., ..]) - &hx.slice(s![..-1, ..])) / dx;
            let b2d = Array2::from_shape_vec((nx-2, 1), b_x_int.clone()).unwrap();
            let a2d = Array2::from_shape_vec((nx-2, 1), a_x_int.clone()).unwrap();
            phi_x.assign(&(&b2d * &phi_x + &a2d * &dh_dx));
            
            &dh_dx + &phi_x
        };

        // === CPML Y Direction (Berechnung der gestreckten Ableitung) ===
        let d2u_dy2_stretched = {
            let du_dy = (&u.slice(s![1..-1, 1..]) - &u.slice(s![1..-1, ..-1])) / dy;
            let b2d = Array2::from_shape_vec((1, ny-1), b_y_half.clone()).unwrap();
            let a2d = Array2::from_shape_vec((1, ny-1), a_y_half.clone()).unwrap();
            psi_y.assign(&(&b2d * &psi_y + &a2d * &du_dy));

            let hy = &du_dy + &psi_y;
            let dh_dy = (&hy.slice(s![.., 1..]) - &hy.slice(s![.., ..-1])) / dy;
            let b2d = Array2::from_shape_vec((1, ny-2), b_y_int.clone()).unwrap();
            let a2d = Array2::from_shape_vec((1, ny-2), a_y_int.clone()).unwrap();
            phi_y.assign(&(&b2d * &phi_y + &a2d * &dh_dy));
            
            &dh_dy + &phi_y
        };

        // === Parallel Wavefield Update (Mit CPML & korrekten Indexgrenzen) ===
        let damp = params.dampfung_global;

        u_next.slice_mut(s![1..-1, 1..-1])
            .axis_iter_mut(Axis(0))
            .into_par_iter()
            .zip(u.slice(s![1..-1, 1..-1]).axis_iter(Axis(0)))
            .zip(u_past.slice(s![1..-1, 1..-1]).axis_iter(Axis(0)))
            .zip(v.slice(s![1..-1, 1..-1]).axis_iter(Axis(0)))
            .zip(d2u_dx2_stretched.axis_iter(Axis(0)))
            .zip(d2u_dy2_stretched.axis_iter(Axis(0)))
            .for_each(|(((((mut next_row, curr_row), past_row), v_row), d2x_row), d2y_row)| {
                // Nun können wir sicher über die gesamte Breite (0..ny-2) iterieren!
                for j in 0..(ny - 2) {
                    let curr = curr_row[j];
                    let past = past_row[j];
                    let vel = v_row[j];
                    let d2x = d2x_row[j];
                    let d2y = d2y_row[j];

                    // Hier fließen jetzt die echten CPML-strikten Ableitungen ein
                    let acc = dt * dt * vel * vel * (d2x + d2y);

                    if damp > 0.0 {
                        let term = (1.0 - 0.5 * damp * dt) * past;
                        next_row[j] = (2.0 * curr - term + acc) / (1.0 + 0.5 * damp * dt);
                    } else {
                        next_row[j] = 2.0 * curr - past + acc;
                    }
                }
            });

        // Source injection
        let t = (i as f64) * dt;
        if t <= 2.5 / params.f_dominant {
            let amp = ricker_wavelet(t, params.f_dominant);
            u_next[[qx, qy]] += amp * dt * dt * v[[qx, qy]] * v[[qx, qy]];
        }

        u_past.assign(&u);
        u.assign(&u_next);

        if i % save_step == 0 {
            write_frame(&mut file, &u, pml)?;
        }
    }

    println!("\n[+] Simulation finished.");
    Ok(())
}

// Fixed memory layout to match row-major order: Y outer loop, X inner loop
fn write_frame(f: &mut File, u: &Array2<f64>, pml: usize) -> io::Result<()> {
    let (nx, ny) = u.dim();
    let mut buf = Vec::with_capacity((nx - 2 * pml) * (ny - 2 * pml));

    for j in pml..ny - pml {       // Outer loop: Y (height)
        for i in pml..nx - pml {   // Inner loop: X (width)
            buf.push(u[[i, j]] as f32);
        }
    }
    f.write_all(bytemuck::cast_slice(&buf))
}
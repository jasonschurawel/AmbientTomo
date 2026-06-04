use crate::Params;
use plotters::prelude::*;
use std::fs::File;
use std::io::{self, Read};

pub fn create_video(bin_path: &str, video_path: &str, _params: &Params) -> io::Result<()> {
    let mut f = File::open(bin_path)?;
    let mut buf = [0u8; 8];

    f.read_exact(&mut buf)?; let _dx = f64::from_le_bytes(buf);
    f.read_exact(&mut buf)?; let _dy = f64::from_le_bytes(buf);
    f.read_exact(&mut buf)?; let _dt = f64::from_le_bytes(buf);
    
    let mut tmp = [0u8; 4];
    f.read_exact(&mut tmp)?; let pml = u32::from_le_bytes(tmp) as usize;
    f.read_exact(&mut tmp)?; let nx = u32::from_le_bytes(tmp) as usize;
    f.read_exact(&mut tmp)?; let ny = u32::from_le_bytes(tmp) as usize;

    let width = nx - 2 * pml;
    let height = ny - 2 * pml;

    // === Physikalisch korrektes Scaling (1 Gitterpunkt = `scale` Pixel in X und Y) ===
    let scale = 1; // Setze auf 2 oder 3, wenn das Video insgesamt höher aufgelöst sein soll
    let out_w = ((width * scale) as u32 / 2) * 2;
    let out_h = ((height * scale) as u32 / 2) * 2;

    println!("Rendering video {}x{} (Aspect Ratio: {:.2}:1) ...", out_w, out_h, (out_w as f32 / out_h as f32));

    let mut frame_data = vec![0f32; width * height];
    let mut frame_count = 0;

    while f.read_exact(bytemuck::cast_slice_mut(&mut frame_data)).is_ok() {
        let frame_filename = format!("frame_{:05}.png", frame_count);
        let backend = BitMapBackend::new(&frame_filename, (out_w, out_h));
        let root = backend.into_drawing_area();

        root.fill(&WHITE).unwrap();

        let max_amp = frame_data.iter().map(|&x| x.abs()).fold(0.0f32, f32::max).max(0.05);

        // Iteration über jeden einzelnen Gitterpunkt ohne Datenverlust
        for y in 0..height {
            for x in 0..width {
                let val = frame_data[y * width + x];
                let intensity = (val / max_amp).clamp(-1.0, 1.0);

                if intensity.abs() > 0.015 {
                    let color = if intensity > 0.0 {
                        RGBColor((intensity * 220.0) as u8, 60, 80)
                    } else {
                        RGBColor(60, 60, (-intensity * 240.0) as u8)
                    };

                    // Koordinaten werden multipliziert statt dividiert -> Kein Clipping mehr!
                    let _ = root.draw(&Rectangle::new(
                        [
                            ((x * scale) as i32, (y * scale) as i32), 
                            (((x + 1) * scale) as i32, ((y + 1) * scale) as i32)
                        ],
                        color.filled(),
                    ));
                }
            }
        }

        root.present().unwrap();
        frame_count += 1;

        if frame_count % 50 == 0 {
            println!("  Rendered {} frames...", frame_count);
        }
    }

    println!("✅ Generated {} frames. Creating video...", frame_count);

    let _ = std::process::Command::new("ffmpeg")
        .args([
            "-y", "-framerate", "60",
            "-i", "frame_%05d.png",
            "-c:v", "libx264", "-pix_fmt", "yuv420p",
            "-preset", "medium", "-crf", "18",
            video_path
        ])
        .status();

    // Cleanup files from disk safely
    for i in 0..frame_count {
        let _ = std::fs::remove_file(format!("frame_{:05}.png", i));
    }

    println!("✅ Video saved: {}", video_path);
    Ok(())
}
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ozean-Simulation (sim1.py)
--------------------------
Simuliert einen 10x3 km Ozean mit einem fraktalen (Perlin-ähnlichen) Meeresboden,
einer perfekten reflektierenden Wasseroberfläche und einem Frequenz-Chirp-Sonar.
"""

import os
import pickle
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.signal import chirp

from wavelib import wellen_simulieren_2d_cpml

if __name__ == "__main__":
    output_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tmp")
    os.makedirs(output_dir, exist_ok=True)
    
    # =========================================================================
    # --- METRISCHE OZEAN-PARAMETER -------------------------------------------
    # =========================================================================
    t_sim_sec = 6.0          # Genug Zeit für Echos vom 3km tiefen Boden
    breite_m = 10000.0       # 10 km breit
    tiefe_m = 4000.0         # 4 km gesamtes Grid (davon ~3km Wasser, ~1km Boden)
    pml_dicke_m = 300.0      # Absorbierende Ränder (außer oben!)
    
    # Materialeigenschaften
    v_wasser = 1500.0        # Schallgeschwindigkeit in Salzwasser (m/s)
    v_fels = 3500.0          # Schallgeschwindigkeit im Basalt/Meeresboden (m/s)
    
    # Dämpfung (Annäherung)
    dampfung_wasser = 0.00   # Wasser dämpft tiefe Frequenzen fast nicht (Verlustfrei)
    dampfung_fels = 0.80     # Fels schluckt enorm viel Energie
    
    # Die U-Boot Sonar-Quelle
    quelle_x_m = 5000.0      # Mittig im Ozean
    quelle_tiefe_m = 500.0   # 500m unter der perfekten Wasseroberfläche
    
    # Das Chirp-Signal (Frequency Sweep von 5 bis 25 Hz)
    f_start = 5.0
    f_end = 25.0
    f_max_effektiv = 30.0    # Für die Gitter-Stabilität (Punkte pro Welle)
    punkte_pro_welle = 8.0   # 8-10 ist ideal für Stabilität vs. RAM
    
    # Rendering
    fps = 60
    zeitlupe_faktor = 3.0    # 6s Physik -> 18s Video

    # =========================================================================
    # --- DYNAMISCHER GITTER-AUFBAU ---
    # =========================================================================
    lambda_min = v_wasser / f_max_effektiv
    d_x = lambda_min / punkte_pro_welle
    d_y = d_x
    
    pml_punkte = int(pml_dicke_m / d_x)
    n_x_innen = int(breite_m / d_x)
    n_y_innen = int(tiefe_m / d_y)
    n_x = n_x_innen + 2 * pml_punkte
    n_y = n_y_innen + 2 * pml_punkte
    
    # 1. Arrays initialisieren
    v = np.ones((n_x, n_y)) * v_wasser
    dampfung_2d = np.ones((n_x, n_y)) * dampfung_wasser
    
    # 2. Fraktalen Meeresboden generieren (Pseudo-Perlin Noise in 1D)
    x_coords = np.linspace(0, breite_m, n_x)
    # Basis-Tiefe: 3000m. Überlagert mit 3 Sinuswellen verschiedener Frequenz/Amplitude
    boden_tiefe_m = 3000.0 - 400.0 * np.sin(x_coords/1200.0) \
                           + 150.0 * np.sin(x_coords/400.0) \
                           -  50.0 * np.cos(x_coords/150.0)
                           
    # 3. Den Boden in die Arrays "stanzen"
    for i in range(n_x):
        boden_idx = int(boden_tiefe_m[i] / d_y) + pml_punkte
        v[i, boden_idx:] = v_fels
        dampfung_2d[i, boden_idx:] = dampfung_fels
        
    # =========================================================================
    # --- SIGNAL-GENERIERUNG ---
    # =========================================================================
    d_t = 0.5 / (np.max(v) * np.sqrt(1.0/d_x**2 + 1.0/d_y**2))
    n_t_total = int(t_sim_sec / d_t)
    
    # Erzeuge ein präzises Chirp-Array für die Engine
    zeit_array = np.linspace(0, t_sim_sec, n_t_total)
    dauer_chirp = 0.5 # Der Sweep dauert 0.5 Sekunden
    chirp_raw = chirp(zeit_array, f0=f_start, t1=dauer_chirp, f1=f_end, method='linear')
    
    # Eine Hanning-Fensterfunktion (Tapering) verhindert numerische Knack-Geräusche
    taper = np.zeros_like(zeit_array)
    taper_idx = int(dauer_chirp / d_t)
    taper[:taper_idx] = np.hanning(taper_idx)
    quell_signal = chirp_raw * taper
    
    # =========================================================================
    # --- SIMULATION STARTEN ---
    # =========================================================================
    video_dauer_sec = t_sim_sec * zeitlupe_faktor
    total_video_frames = int(video_dauer_sec * fps)
    speicher_schritt = max(1, int(n_t_total / total_video_frames))
    
    q_x = int(quelle_x_m / d_x) + pml_punkte
    q_y = int(quelle_tiefe_m / d_y) + pml_punkte

    def finde_naechste_nummer(ordner):
        Dateien = [f for f in os.listdir(ordner) if f.startswith("sim_") and f.endswith(".mp4")]
        if not Dateien: return 1
        nummern = [int(f.split("_")[1]) for f in Dateien if len(f.split("_")) > 1 and f.split("_")[1].isdigit()]
        return max(nummern) + 1 if nummern else 1
    
    naechste_id = finde_naechste_nummer(output_dir)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    pkl_filename = os.path.join(output_dir, f"sim_{naechste_id:04d}_{timestamp}.pkl")
    video_filename = os.path.join(output_dir, f"sim_{naechste_id:04d}_{timestamp}.mp4")
    
    print(f"\n[*] Ozean-Dimensionen: {breite_m/1000} km breit, {tiefe_m/1000} km tief.")
    print(f"[*] Grid: {n_x}x{n_y} Punkte (dx={d_x:.2f}m). RAM-Streaming aktiv.\n")
    
    # --- PHASE 1: RECHNEN ---
    wellen_simulieren_2d_cpml(v, dampfung_2d, quell_signal, t_sim_sec, d_t, d_x, d_y, 
                              pml_punkte, q_x, q_y, f_max_effektiv, speicher_schritt, pkl_filename)
    
    # --- PHASE 2: RENDERING ---
    print(f"\n[+] Erzeuge Visualisierung: Sonar im tiefen Ozean...")
    fig, ax = plt.subplots(figsize=(10, 5))
    
    with open(pkl_filename, 'rb') as f_pkl:
        # --- DYNAMISCHES RENDERING ---
        # Berechne das physikalische Verhältnis für die Figure-Größe
        verhaeltnis = breite_m / tiefe_m
        fig, ax = plt.subplots(figsize=(10, 10 / verhaeltnis))
        
        # Header ignorieren
        pickle.load(f_pkl) 
        erster_frame = pickle.load(f_pkl)
        u_echt = erster_frame[pml_punkte:-pml_punkte, pml_punkte:-pml_punkte]
        
        im = ax.imshow(u_echt.T, cmap='seismic', vmin=-0.05, vmax=0.05, 
               extent=[0, breite_m/1000, tiefe_m/1000, 0], 
               aspect='equal')
        
        # Zeichne den fraktalen Meeresboden nach
        ax.plot(x_coords/1000, boden_tiefe_m/1000, color='black', linewidth=2, label="Meeresboden")
        ax.plot([quelle_x_m/1000], [quelle_tiefe_m/1000], 'y^', markersize=10, label="U-Boot (Sonar)")
        
        # Die Wasseroberfläche
        ax.axhline(0, color='cyan', linestyle='-', linewidth=3, label="Wasseroberfläche")
        
        ax.set_xlabel("Distanz (km)")
        ax.set_ylabel("Tiefe (km)")
        ax.legend(loc="upper right", facecolor='white', framealpha=0.8)
        fig.colorbar(im, label="Akkustischer Druck", pad=0.02)
        fig.tight_layout()
        
        writer = animation.FFMpegWriter(fps=fps, bitrate=5000)
        
        frame_zaehler = 0
        with writer.saving(fig, video_filename, dpi=120):
            try:
                while True:
                    naechster_frame = pickle.load(f_pkl)
                    u_echt = naechster_frame[pml_punkte:-pml_punkte, pml_punkte:-pml_punkte]
                    
                    # Dynamischer Kontrast
                    std_dev = np.std(u_echt)
                    grenzwert = std_dev * 4.0 if std_dev > 1e-6 else 0.001
                    
                    im.set_array(u_echt.T)
                    im.set_clim(vmin=-grenzwert, vmax=grenzwert)
                    
                    aktuelle_zeit_s = frame_zaehler * speicher_schritt * d_t
                    ax.set_title(f"10km Ozean Sonar (Chirp {f_start}-{f_end} Hz) | Zeit: {aktuelle_zeit_s:.2f} s")
                    
                    writer.grab_frame()
                    frame_zaehler += 1
                    
                    if frame_zaehler % 50 == 0:
                        print(f"\r    Rendering: {frame_zaehler} Frames geschrieben.", end="", flush=True)
            except EOFError:
                print(f"\r    Rendering: [100%] Abgeschlossen. Video gespeichert unter {video_filename}")

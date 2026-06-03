#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Seismische 2D-Wellensimulation - Bazel Execution Script (simulation.py)
---------------------------------------------------------------------
Dieses Skript ist speziell für den Aufruf über das Bazel Build-System optimiert.
Es trennt Physik und Rendering über eine transiente PKL-Datei und schreibt
das finale Video direkt an den von Bazel vorgegebenen Zielort.
"""

import os
import sys
import pickle
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from wellen_simulation import wellen_simulieren_2d_cpml

if __name__ == "__main__":
    # Bazel übergibt den gewünschten Pfad für das Output-Video als erstes Argument ($@)
    if len(sys.argv) < 2:
        print("Fehler: Dieses Skript benötigt den Ausgabepfad für das Video als Argument.")
        print("Verwendung: python3 simulation.py <pfad_zum_video.mp4>")
        sys.exit(1)
        
    video_filename = sys.argv[1]
    # Die temporäre PKL-Datei legen wir parallel zur Ausgabedatei ab
    pkl_filename = video_filename.replace(".mp4", ".pkl")

    # =========================================================================
    # --- METRISCHE PHYSIK-PARAMETER (In Metern & Sekunden) ---
    # =========================================================================
    t_sim_sec = 2.5          # Physikalische Simulationsdauer (Sekunden)
    breite_m = 5000.0        # Horizontale Ausdehnung des Kerngebiets (Meter)
    tiefe_m = 1000.0         # Vertikale Ausdehnung des Kerngebiets (Meter)
    grenzschicht_m = 500.0   # Tiefe der geologischen Schichtgrenze (Meter)
    pml_dicke_m = 250.0      # Dicke des absorbierenden CPML-Rands (Meter)
    
    quelle_x_m = 400.0       # Position der seismischen Quelle auf der x-Achse (Meter)
    quelle_tiefe_m = 100.0   # Tiefe der Quelle (Meter)
    
    v_oben = 1500.0          # Schallgeschwindigkeit obere Schicht (m/s)
    v_unten = 3000.0         # Schallgeschwindigkeit untere Schicht (m/s)
    
    f_dominant = 15.0        # Dominante Grundfrequenz der Quelle (Hz)
    f_max_effektiv = 40.0    # Eingestellte Oberton-Grenze für die Auflösung (Hz)
    punkte_pro_welle = 10.0  # Auflösung: Gitterpunkte pro kleinster Wellenlänge
    dampfung_global = 0.4    # Reibungskoeffizient für natürliche Mediendämpfung

    # =========================================================================
    # --- VISUELLE RENDERING-PARAMETER (Entkoppelt) -----------
    # =========================================================================
    fps = 60                 # Bilder pro Sekunde im Video
    zeitlupe_faktor = 6.0    # Zeitdehnung (2.5s Physik -> 15s Video)
    
    video_dauer_sec = t_sim_sec * zeitlupe_faktor
    total_video_frames = int(video_dauer_sec * fps)

    # =========================================================================
    # --- DYNAMISCHE STRUKTUR- UND STABILITÄTSANALYSE ---------
    # =========================================================================
    lambda_min = v_oben / f_max_effektiv
    d_x = lambda_min / punkte_pro_welle
    d_y = d_x
    
    pml_punkte = int(pml_dicke_m / d_x)
    n_x_innen = int(breite_m / d_x)
    n_y_innen = int(tiefe_m / d_y)
    n_x = n_x_innen + 2 * pml_punkte
    n_y = n_y_innen + 2 * pml_punkte
    
    v = np.ones((n_x, n_y)) * v_oben
    grenzschicht_idx = int(grenzschicht_m / d_y) + pml_punkte
    v[:, grenzschicht_idx:] = v_unten
    
    d_t = 0.5 / (np.max(v) * np.sqrt(1.0/d_x**2 + 1.0/d_y**2))
    n_t_total = int(t_sim_sec / d_t)
    
    speicher_schritt = max(1, int(n_t_total / total_video_frames))
    
    q_x = int(quelle_x_m / d_x) + pml_punkte
    q_y = int(quelle_tiefe_m / d_y) + pml_punkte
    
    print(f"\n=========================================================================")
    print(f"BAZEL MONOREPO GEOPHYSIK PIPELINE")
    print(f"=========================================================================")
    print(f"[*] Gitterweite (dx, dy):    {d_x:.2f} m  (Kritische Wellenlänge: {lambda_min:.2f} m)")
    print(f"[*] Interner Zeitschritt:    {d_t:.6f} s  (Gesamtschritte: {n_t_total})")
    print(f"[*] Modell-Dimensionen:      {n_x} x {n_y} Gitterpunkte (inkl. {pml_punkte} PML-Zellen)")
    print(f"[*] Video-Spezifikation:     {video_dauer_sec:.1f}s Dauer bei {fps} FPS")
    print(f"[*] Zielpfad (Bazel Sandbox): {video_filename}")
    print(f"=========================================================================\n")
    
    # --- PHASE 1 - SIMULATION (DIRECT-TO-DISK) ---
    wellen_simulieren_2d_cpml(v, t_sim_sec, d_t, d_x, d_y, pml_punkte, q_x, q_y, 
                              f_dominant, dampfung_global, speicher_schritt, pkl_filename)
    
    # --- PHASE 2 - RENDERING (STREAMING AUS PKL) ---
    print(f"\n[+] Starte Phase 2: Videorendering über sequentiellen PKL-Import...")
    fig, ax = plt.subplots(figsize=(10, 4.5))
    
    with open(pkl_filename, 'rb') as f_pkl:
        header = pickle.load(f_pkl)
        erster_frame = pickle.load(f_pkl)
        u_echt = erster_frame[pml_punkte:-pml_punkte, pml_punkte:-pml_punkte]
        
        im = ax.imshow(u_echt.T, cmap='bwr', vmin=-0.05, vmax=0.05, 
                       extent=[0, breite_m/1000, tiefe_m/1000, 0])
        
        ax.axhline(grenzschicht_m/1000, color='black', linestyle='--', alpha=0.6, label="Grenzschicht (500m)")
        ax.plot([quelle_x_m/1000], [quelle_tiefe_m/1000], 'y*', markersize=11, label="Seismische Quelle")
        ax.set_xlabel("Horizontale Distanz (km)")
        ax.set_ylabel("Tiefe (km)")
        ax.legend(loc="lower right")
        fig.colorbar(im, label="Amplitude (Dynamische AGC)", pad=0.02)
        fig.tight_layout()
        
        writer = animation.FFMpegWriter(fps=fps, bitrate=4500)
        frame_zaehler = 0
        
        with writer.saving(fig, video_filename, dpi=110):
            writer.grab_frame()
            frame_zaehler += 1
            
            try:
                while True:
                    naechster_frame = pickle.load(f_pkl)
                    u_echt = naechster_frame[pml_punkte:-pml_punkte, pml_punkte:-pml_punkte]
                    
                    std_dev = np.std(u_echt)
                    grenzwert = std_dev * 3.0 if std_dev > 1e-5 else 0.002
                    
                    im.set_array(u_echt.T)
                    im.set_clim(vmin=-grenzwert, vmax=grenzwert)
                    
                    aktuelle_zeit_s = frame_zaehler * speicher_schritt * d_t
                    ax.set_title(f"Seismische Wellenfront - Zeit: {aktuelle_zeit_s:.3f} s")
                    
                    writer.grab_frame()
                    frame_zaehler += 1
                    
                    if frame_zaehler % max(1, total_video_frames // 20) == 0:
                        print(f"\r    Fortschritt Rendering: {frame_zaehler} Frames verarbeitet.", end="", flush=True)
                        
            except EOFError:
                print(f"\r    Fortschritt Rendering: [100%] Alle {frame_zaehler} Frames verarbeitet.")
                
    print(f"[+] Video erfolgreich exportiert: {video_filename}")
    plt.close()
    
    # --- SAUBERES BAZEL CONTEXT CLEANUP ---
    if os.path.exists(pkl_filename):
        os.remove(pkl_filename)
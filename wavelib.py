#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Wellen-Simulation Backend (wellen_simulation.py)
------------------------------------------------
Dieses Modul stellt die Kernfunktionalität zur Simulation akustischer Wellen 
in einem 2D-Medium bereit. Es nutzt die Finite-Differenzen-Methode (FDTD) 
auf einem versetzten Gitter (Staggered Grid) und implementiert hochwirksame 
absorbierende Randbedingungen (Convolutional Perfectly Matched Layer - CPML).

Optimierung:
Die Ergebnisse werden Frame für Frame direkt in eine Pickle-Datei (.pkl) 
gestreamt. Dadurch verbleibt zu keinem Zeitpunkt das gesamte 3D-Wellenfeld 
im Arbeitsspeicher (RAM-Schutz).
"""

import numpy as np
import pickle

def erzeuge_cpml_koeffizienten(n, pml_punkte, d, v_max, d_t, f_dominant):
    """
    Berechnet die CPML-Koeffizienten (a und b) für eine Dimension.
    Es wird strikt zwischen ganzzahligen Knoten (int) und Halbschritten (half)
    des Staggered Grids unterschieden, um künstliche Reflexionen zu minimieren.
    """
    R_theoretisch = 1e-6
    L = pml_punkte * d
    
    sigma_max = -(3.0 * v_max * np.log(R_theoretisch)) / (2.0 * L)
    alpha_max = np.pi * f_dominant
    
    sigma_int = np.zeros(n - 2)
    alpha_int = np.zeros(n - 2)
    
    for idx in range(n - 2):
        g_idx = idx + 1  
        if g_idx < pml_punkte:
            dist = (pml_punkte - g_idx) * d
            sigma_int[idx] = sigma_max * (dist / L)**2
            alpha_int[idx] = alpha_max * (1.0 - dist / L)
        elif g_idx >= n - pml_punkte:
            dist = (g_idx - (n - 1 - pml_punkte)) * d
            sigma_int[idx] = sigma_max * (dist / L)**2
            alpha_int[idx] = alpha_max * (1.0 - dist / L)
            
    sigma_half = np.zeros(n - 1)
    alpha_half = np.zeros(n - 1)
    
    for idx in range(n - 1):
        h_idx = idx + 0.5  
        if h_idx < pml_punkte:
            dist = (pml_punkte - h_idx) * d
            sigma_half[idx] = sigma_max * (dist / L)**2
            alpha_half[idx] = alpha_max * (1.0 - dist / L)
        elif h_idx > n - 1 - pml_punkte:
            dist = (h_idx - (n - 1 - pml_punkte)) * d
            sigma_half[idx] = sigma_max * (dist / L)**2
            alpha_half[idx] = alpha_max * (1.0 - dist / L)
            
    b_int = np.exp(-(sigma_int + alpha_int) * d_t)
    a_int = np.zeros_like(sigma_int)
    mask_int = (sigma_int + alpha_int) > 0
    a_int[mask_int] = (b_int[mask_int] - 1.0) * sigma_int[mask_int] / (sigma_int[mask_int] + alpha_int[mask_int])
    
    b_half = np.exp(-(sigma_half + alpha_half) * d_t)
    a_half = np.zeros_like(sigma_half)
    mask_half = (sigma_half + alpha_half) > 0
    a_half[mask_half] = (b_half[mask_half] - 1.0) * sigma_half[mask_half] / (sigma_half[mask_half] + alpha_half[mask_half])
    
    return a_int, b_int, a_half, b_half

def ricker_wavelet(t, f_dominant):
    """ Generiert die Quell-Amplitude (Druckimpuls) zum Zeitpunkt t. """
    t0 = 1.2 / f_dominant  
    tau = np.pi * f_dominant * (t - t0)
    return (1.0 - 2.0 * tau**2) * np.exp(-tau**2)

def wellen_simulieren_2d_cpml(v, t_sim_sec, d_t, d_x, d_y, pml_punkte, q_x, q_y, f_dominant, dampfung_global, speicher_schritt, pkl_filepath):
    """
    Führt die 2D-Wellensimulation durch und streamt die Daten direkt auf die Festplatte.
    """
    n_x, n_y = v.shape
    n_t_total = int(t_sim_sec / d_t)
    
    a_x_int, b_x_int, a_x_half, b_x_half = erzeuge_cpml_koeffizienten(n_x, pml_punkte, d_x, np.max(v), d_t, f_dominant)
    a_y_int, b_y_int, a_y_half, b_y_half = erzeuge_cpml_koeffizienten(n_y, pml_punkte, d_y, np.max(v), d_t, f_dominant)
    
    a_x_int, b_x_int = a_x_int[:, np.newaxis], b_x_int[:, np.newaxis]
    a_x_half, b_x_half = a_x_half[:, np.newaxis], b_x_half[:, np.newaxis]
    a_y_int, b_y_int = a_y_int[np.newaxis, :], b_y_int[np.newaxis, :]
    a_y_half, b_y_half = a_y_half[np.newaxis, :], b_y_half[np.newaxis, :]
    
    # RAM-Schutz: Es existieren immer nur drei 2D-Zeitschichten gleichzeitig im Speicher
    u = np.zeros((n_x, n_y))
    u_past = np.zeros((n_x, n_y))
    u_next = np.zeros((n_x, n_y))
    
    psi_x = np.zeros((n_x - 1, n_y - 2))
    phi_x = np.zeros((n_x - 2, n_y - 2))
    psi_y = np.zeros((n_x - 2, n_y - 1))
    phi_y = np.zeros((n_x - 2, n_y - 2))
    
    # Datei öffnen und ersten Header schreiben
    with open(pkl_filepath, 'wb') as f_pkl:
        metadata = {
            'v': v, 'd_x': d_x, 'd_y': d_y, 'd_t': d_t,
            'pml_punkte': pml_punkte, 'speicher_schritt': speicher_schritt,
            't_sim_sec': t_sim_sec, 'n_t_total': n_t_total
        }
        pickle.dump(metadata, f_pkl)
        
        print(f"[+] Starte Physik-Engine ({n_t_total} Schritte). Stream läuft in PKL-Datei...")
        
        frames_saved = 0
        for i in range(n_t_total - 2):
            if i % max(1, n_t_total // 20) == 0 or i == n_t_total - 3:
                prozent = int((i + 2) / n_t_total * 100)
                print(f"\r    Simuliere: [{prozent:03d}%] Schritt {i+2}/{n_t_total}", end="", flush=True)
                
            # CPML-Berechnungen ...
            du_dx = (u[1:, 1:-1] - u[:-1, 1:-1]) / d_x
            psi_x = b_x_half * psi_x + a_x_half * du_dx
            Hx = du_dx + psi_x
            dH_dx = (Hx[1:, :] - Hx[:-1, :]) / d_x
            phi_x = b_x_int * phi_x + a_x_int * dH_dx
            d2u_dx2_stretched = dH_dx + phi_x
            
            du_dy = (u[1:-1, 1:] - u[1:-1, :-1]) / d_y
            psi_y = b_y_half * psi_y + a_y_half * du_dy
            Hy = du_dy + psi_y
            dH_dy = (Hy[:, 1:] - Hy[:, :-1]) / d_y
            phi_y = b_y_int * phi_y + a_y_int * dH_dy
            d2u_dy2_stretched = dH_dy + phi_y
            
            beschleunigung = (d_t**2) * (v[1:-1, 1:-1]**2) * (d2u_dx2_stretched + d2u_dy2_stretched)
            
            if dampfung_global > 0.0:
                term_past = (1.0 - 0.5 * dampfung_global * d_t) * u_past[1:-1, 1:-1]
                u_next[1:-1, 1:-1] = (2.0 * u[1:-1, 1:-1] - term_past + beschleunigung) / (1.0 + 0.5 * dampfung_global * d_t)
            else:
                u_next[1:-1, 1:-1] = 2.0 * u[1:-1, 1:-1] - u_past[1:-1, 1:-1] + beschleunigung
            
            zeit = i * d_t
            if zeit <= (2.5 / f_dominant):
                amplituden_wert = ricker_wavelet(zeit, f_dominant)
                u_next[q_x, q_y] += amplituden_wert * (d_t**2) * (v[q_x, q_y]**2)
            
            u_past[:, :] = u
            u[:, :] = u_next
            
            # DIRECT-TO-DISK STREAMING (RAM-SCHUTZ)
            if i % speicher_schritt == 0:
                pickle.dump(u_next.astype(np.float32), f_pkl)
                frames_saved += 1
                
    print(f"\n[+] Physik beendet. {frames_saved} Frames erfolgreich in PKL gestreamt.")
# AmbientTomo – Synthetic Aperture Sonar & Tomography

The primary goal of this project is the development of an imaging Synthetic Aperture Sonar (SAS) system utilizing modern tomographic reconstruction techniques. Currently, development is focused on the implementation and validation of the Forward Modeling process.

## What is happening here?

* **Forward Modeling:** To test a sonar system or tomographic algorithms, physical reality must be accurately simulated in a computer environment. To achieve this, the system calculates the propagation of acoustic waves (initially in a two-dimensional medium). An artificial signal (pulse) excites the medium, propagates through it, refracts at material boundaries, and generates reflections. These synthetic data points form the essential foundation for subsequent inversion back into a high-resolution image.
* **Absorbing Boundary Conditions:** To accurately represent infinite domains (such as the open ocean or deep geological structures) within a limited digital grid, specialized boundaries are used to capture wavefronts. Artificial echoes at the edges of the simulation area are suppressed, ensuring that signals exit the frame cleanly without distorting the sensor data. The simulation utilizes the Finite-Difference Time-Domain (FDTD) method on a staggered grid and implements highly effective absorbing boundary conditions (Convolutional Perfectly Matched Layer - CPML).
* **Visualization & Dynamic Gain Control:** Since acoustic amplitudes decay significantly with increasing distance from the signal source, an automatic gain control mechanism ensures that both the high-energy initial pulse and extremely faint deep reflections are rendered clearly and with high contrast in the visual output.

## Integration with Bazel
The entire system is built to be fully deterministic. Any modification to the physical forward modeling or the wave propagation parameters triggers an automatically reproducible recalculation of the simulation results. Unchanged states are loaded directly from the global build cache without requiring additional computation.

(Of course it it possible to run the simulation without bazel by executing the simulation.py)
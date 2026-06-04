# AmbientTomo – Synthetic Aperture Sonar & Tomography

The primary goal of this project is the development of an imaging Synthetic Aperture Sonar (SAS) system utilizing modern tomographic reconstruction techniques. Currently, development is focused on the implementation and validation of the Forward Modeling process.

## What is happening here?

* **Forward Modeling:** To test a sonar system or tomographic algorithms, physical reality must be accurately simulated in a computer environment. To achieve this, the system calculates the propagation of acoustic waves (initially in a two-dimensional medium). An artificial signal (pulse) excites the medium, propagates through it, refracts at material boundaries, and generates reflections. These synthetic data points form the essential foundation for subsequent inversion back into a high-resolution image.
* **Absorbing Boundary Conditions:** To accurately represent infinite domains (such as the open ocean or deep geological structures) within a limited digital grid, specialized boundaries are used to capture wavefronts. Artificial echoes at the edges of the simulation area are suppressed, ensuring that signals exit the frame cleanly without distorting the sensor data. The simulation utilizes the Finite-Difference Time-Domain (FDTD) method on a staggered grid and implements highly effective absorbing boundary conditions (Convolutional Perfectly Matched Layer - CPML).
* **Visualization & Dynamic Gain Control:** Since acoustic amplitudes decay significantly with increasing distance from the signal source, an automatic gain control mechanism ensures that both the high-energy initial pulse and extremely faint deep reflections are rendered clearly and with high contrast in the visual output.

---

## Governing Physical Equations

The forward model solves the **2D Acoustic Wave Equation** with a global damping term to model natural medium attenuation:

$$\frac{\partial^2 u}{\partial t^2} + \gamma \frac{\partial u}{\partial t} = v(x,y)^2 \left( \frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} \right) + s(x,y,t)$$

Where:
* $u(x,y,t)$ is the wavefield amplitude (pressure/displacement).
* $v(x,y)$ is the spatial medium velocity (e.g., $1500\text{ m/s}$ top layer, $3000\text{ m/s}$ bottom layer).
* $\gamma$ is the global damping coefficient (`dampfung_global`).
* $s(x,y,t)$ is the seismic/sonar source term (injected as a Ricker Wavelet).

---

## Future Roadmap: Higher Spatial Stencil Orders (FDTD 4th Order)

To drastically increase numerical accuracy while optimizing computing time, a planned feature is upgrading the spatial finite-difference operators from the current 2nd-order scheme to a **4th-order** accurate stencil.

### 1. Mathematical Finite-Difference Approximations
While the current 2nd-order scheme uses a 3-point central difference, a 4th-order spatial stencil expands the operator to 5 points to approximate the second-order partial spatial derivatives:

$$\frac{\partial^2 u}{\partial x^2} \approx \frac{-u_{i+2,j} + 16u_{i+1,j} - 30u_{i,j} + 16u_{i-1,j} - u_{i-2,j}}{12 \cdot \Delta x^2}$$

$$\frac{\partial^2 u}{\partial y^2} \approx \frac{-u_{i,j+2} + 16u_{i,j+1} - 30u_{i,j} + 16u_{i,j-1} - u_{i,j-2}}{12 \cdot \Delta y^2}$$

### 2. Benefits of the 4th-Order Upgrade
* **Drastic Reduction in Numerical Dispersion:** In 2nd-order grids, high-frequency wave components suffer from artificial phase velocity shifts (short wavelengths travel slower than physically correct). A 4th-order stencil suppresses this grid-induced numerical dispersion almost entirely.
* **Coarser Grids ($\Delta x, \Delta y$):** A 2nd-order scheme typically requires $10$ to $20$ grid points per minimum wavelength ($\lambda_{\text{min}}$) to remain accurate. With 4th-order accuracy, this constraint drops to just $4$ to $6$ points per wavelength.
* **Massive Performance Gain:** Increasing $\Delta x$ means the total grid size ($N_x \times N_y$) scales down quadratically. Even though each grid point requires slightly more arithmetic operations, the overall simulation runs multiple times faster and consumes significantly less RAM.

### 3. CPML Considerations
Upgrading the core wave equation requires a synchronized change in the CPML boundaries. Because CPML relies on *first-order* spatial derivatives to advance its memory variables ($\psi$ and $\phi$), using mixed orders would introduce artificial impedance mismatches and reflections at the PML interfaces.

* **4th-Order CPML Derivatives:** First-order spatial derivatives inside the PML regions must be upgraded to a matching 4th-order central difference:
  $$\frac{\partial u}{\partial x} \approx \frac{-u_{i+2,j} + 8u_{i+1,j} - 8u_{i-1,j} + u_{i-2,j}}{12 \cdot \Delta x}$$
* **Staggered Grid Interlocking:** The coefficients must be meticulously re-mapped to account for the half-step grid offsets ($\pm 1/2, \pm 3/2$) typical for staggered-grid FDTD layouts to ensure perfectly non-reflective boundaries.

---

## Rust Implementation (`rust_src`)

The Rust port is designed for maximum throughput using `ndarray` and `rayon` multi-threading. Transitioning to a 4th-order stencil introduces specific implementation changes:

* **Ghost/Boundary Cell Adjustments:** The current 2nd-order inner loop slices the grid at `s![1..-1, 1..-1]`. For a 4th-order stencil, the active compute domain must shrink to `s![2..-2, 2..-2]` to avoid out-of-bounds errors when accessing $i \pm 2$ and $j \pm 2$. Consequently, the outer PML boundary padding might need to be increased slightly.
* **Cache Locality & Multi-Threading:** The parallelized row update loop (`.axis_iter_mut(Axis(0)).into_par_iter()`) remains highly efficient with `rayon`. However, since a 4th-order stencil reads from two adjacent rows up ($i-2$) and two rows down ($i+2$), cache line utilization must be closely monitored to prevent cache thrashing during parallel processing.
* **SIMD Auto-Vectorization:** Rust’s compiler (`rustc` via LLVM) excels at auto-vectorizing fixed stencil weights when compiled with target CPU features flags (e.g., `-C target-cpu=native`). The explicit unrolling of the 5-point stencil coefficients (`-1, 16, -30, 16, -1`) will map cleanly to AVX2/AVX-512 vector hardware instructions, maximizing instructions-per-cycle (IPC).

---

## Integration with Bazel

The entire system is built to be fully deterministic. Any modification to the physical forward modeling or the wave propagation parameters triggers an automatically reproducible recalculation of the simulation results. Unchanged states are loaded directly from the global build cache without requiring additional computation.

*(Of course, it is possible to run the simulation without Bazel by executing `simulation.py` directly).*
FMX — Fast Free‑Molecular Aerodynamics (100–600 km)

Overview
- Panel‑based C++20 library and CLI to compute aerodynamic forces and torques on satellites in LEO (free‑molecular to transition regime).
- Gas–surface interaction: Sentman closed‑form with energy accommodation; optional numerical quadrature for verification.
- Geometry: OBJ/ASCII STL loader, BVH occlusion, parallel per‑facet integration (OpenMP).
- Atmosphere: NRLMSIS2.1 (T, species) + HWM14 (winds) wrappers; stub fallback.
- Validation: unit tests for Sentman behavior, occlusion, angle trends, torque; harness for NASA CR‑313 reference cases.

Build
- Prereqs: CMake ≥ 3.16, C++20 compiler; optional Fortran (gfortran) for local models; OpenMP optional.
- Configure (with local NRLMSIS2.1 + HWM14):
  cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DFMX_ENABLE_OPENMP=ON -DFMX_USE_LOCAL_MODELS=ON
  cmake --build build -j
- Run tests:
  ctest --test-dir build --output-on-failure

CLI
- Help: ./build/fmx_cli --help
- Validation presets:
  - Plate: ./build/fmx_cli --validate plate
  - Two plates (occlusion): ./build/fmx_cli --validate two-plates
  - Cube: ./build/fmx_cli --validate cube
  - Torque (plate, offset CG): ./build/fmx_cli --validate torque-plate
- From JSON config:
  ./build/fmx_cli --config examples/config.sample.json --out out.json
- Output JSON includes F, M, and diagnostics (facets, shadowed, T, tau, Ma range).

Config Schema (minimal)
- geometry: path to mesh (OBJ or ASCII STL)
- cg: [x,y,z] center of gravity (m)
- materials.default: { alpha_E, Tw_K }
- atmosphere:
  - model: "Stub" | "NRLMSIS2" | "NRLMSIS2+HWM14"
  - indices: { F10_7, F10_7A, Kp, Ap (number or array[7]), Ap_daily, Ap_now }
- state: { alt_km, lat_deg, lon_deg, utc (ISO8601 Z), V_sat_mps: [vx,vy,vz] }

Atmosphere Models
- Local models (Fortran) from atm/models/:
  - NRLMSIS2.1 (jacobwilliams/NRLMSIS2.1) — needs msis21.parm at runtime (auto‑copied).
  - HWM14 (gemini3d/hwm14) — needs hwm123114.bin, dwm07b104i.dat, gd2qd.dat (auto‑copied).
- Wrappers provide C bindings (atm/wrappers/*.f90) and C++ interfaces:
  - NRLMSIS2Atmosphere — T and species (mass density via number density × species mass)
  - HWM14 — neutral winds
  - CombinedAtmosphere — NRLMSIS2 + HWM14

Gas–Surface Interaction
- Sentman closed‑form traction coefficients C_N, C_T vs angle, speed ratio; numerically stable (erfc/exp) branches.
- Reflected temperature via energy accommodation (Tuttas et al., 2025): tau = (1-α_E)·E_i/(2kT_i) + α_E·(T_w/T_i).
- Numerical quadrature (Gauss–Hermite) retained as verification path.

Occlusion & Solver
- BVH occluder (median split) with slab AABB and Möller–Trumbore any‑hit.
- Per‑facet parallel integration (OpenMP) with reductions; optional serial path.

Validation Suite
- Unit tests (ctest):
  - sentman_basic — edge cases (normal/grazing; Ma behavior)
  - occlusion_basic — rear plate shadowing
  - validation_trends — angle of attack trends
  - torque_plate_offset — torque consistency Mz ≈ r×F
  - cube_symmetry_drag — symmetry and drag direction checks
  - cr313_reference — harness to compare against NASA CR‑313 reference cases
- CR‑313:
  - Place a file at examples/cr313/cases.json describing cases with expected F/M and tolerance.
  - The test harness will run available cases and compare; if file absent, it skips.

Notes & Limits
- Transition regime blending (R13) is stubbed for now; free‑molecular core is implemented.
- HWM14 and MSIS require data files present in CWD; the wrappers copy them on first use from atm/models/.
- UTC parsing expects ISO8601 Z; timezone offsets are not parsed.

License
- See LICENSE files of third‑party models under atm/models/; FMX project license TBD by the repository owner.


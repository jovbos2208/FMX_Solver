# **FMX: Ein schneller, paralleler Standard für aerodynamische Kräfte & Momente im LEO (100–600 km)**  
*Entwurf eines C++-Projekts (mit Parallelisierung) samt Theorie, Gleichungen, Validierung & Quellen. Diese Datei ist dafür gedacht, mit einem Code-Assistenten Schritt für Schritt eine performante Bibliothek/CLI aufzubauen.*

---

## Ziel & Scope

Wir bauen eine **deterministische, panelbasierte** Methode zur Berechnung von **Kräften und Momenten** durch die Atmosphäre auf beliebigen Satellitengeometrien im Bereich **100–600 km**. Die Methode soll **schneller als DSMC** und im **frei-molekularen bis Übergangsregime** hinreichend genau sein. Kernideen:

- **Atmosphären-State** aus Empirikmodellen (NRLMSIS/DTM/JB2008) + Winde (HWM14).  
- **Gas-Oberflächen-Wechselwirkung (GSI)** via Sentman/CLL-Familie (Cercignani–Lampis–Lord) mit vorab generierten **Koeffizienten-Kernen** (Tabellen/NN-Surrogate).  
- **Panel-Ray-Tracing** zur Okklusion (Schatten) + **massive Parallelisierung** (OpenMP/TBB; optional Embree/CUDA).  
- **Regime-Adapter** für 100–150 km über **Momentsysteme (R13)** als korrigierender Blend-Term.

---

## 1) Physikalischer Hintergrund

### 1.1 Regime (LEO 100–600 km)
- **Knudsen-Zahl** groß → ab ~150–200 km weitgehend **frei-molekular**, 100–150 km **Übergang**.  
- **Zusammensetzung**: ~100–150 km gemischt (N₂/O₂/O), 150–400 km **O** dominant, >400 km zunehmend **He/H**.  
- **Variabilität**: Dichte ~10⁻⁴–10⁻¹² kg m⁻³; starke Abhängigkeit von **F10.7** (Solarflux) & **Kp/Ap** (geomagnetisch).

### 1.2 Atmosphärenmodelle (für die Pipeline)
- **NRLMSIS 2.0** (T, nₛ, ρ), **DTM2013/DTM2020** (T, ρ, Zusammensetzung), **JB2008** (ρ); **HWM14** (Wind).  
- Hinweise zur Güte/UQ: DTM2020 liefert UQ & reduzierte Biases ggü. Vorgängern.

### 1.3 Gas-Oberflächen-Wechselwirkung (GSI)
- **Sentman (1961)**: klassische FM-Koeffizienten (normal/tangential) mit Energie-/Impuls-Akk. und Wandtemperatur.  
- **Cercignani–Lampis (CL/CLL)**: realistischere Streukerne (getrennte Normal-/Tangential-Akk.). Erweiterungen & Vergleiche zu Simulationen.

### 1.4 Übergangsregime-Korrekturen
- **Regularisierte 13-Momenten (R13)** als schlanker Korrektur-Ansatz gegenüber Boltzmann/DSMC, stabiler als Burnett.

---

## 2) Mathematisches Modell

**Relativgeschwindigkeit** \(\mathbf{c}=\mathbf{V}_\text{sat}-\mathbf{u}_\text{wind}\).

Für jede Facette \(i\) (Fläche \(A_i\), Normalen \(\mathbf{n}_i\)) und Spezies \(s\) (Dichte \(\rho_s\), Masse \(m_s\)):

\[
p_{\infty,s}=\rho_s\|\mathbf{c}\|^2,\quad
\theta_i=\arccos\!\big(-\hat{\mathbf{c}}\!\cdot\!\mathbf{n}_i\big),\quad
\mathrm{Ma}_s=\frac{\|\mathbf{c}\|}{\sqrt{kT/m_s}},\quad
\tau=\frac{T_w}{T}
\]

\[
\mathrm{d}\mathbf{F}_{i,s}=p_{\infty,s}A_i\!\left[
C_{N,s}(\theta_i,\mathrm{Ma}_s,\tau,\alpha)\,\mathbf{n}_i
+ C_{T,s}(\theta_i,\mathrm{Ma}_s,\tau,\alpha)\,\hat{\mathbf{t}}_i
\right]
\]

\(C_{N,s},C_{T,s}\) stammen aus **GSI-Kernen** (Sentman/CLL-basierte Integrale → offline tabelliert/NN-Surrogat). Gesamtkraft/-moment:

\[
\mathbf{F}=\sum_{i,s}\mathrm{d}\mathbf{F}_{i,s},\qquad
\mathbf{M}=\sum_i(\mathbf{r}_i-\mathbf{r}_\mathrm{CG})\times\sum_s\mathrm{d}\mathbf{F}_{i,s}
\]

**Okklusion**: Beitrag 0 falls Sichtstrahl entlang \(-\hat{\mathbf{c}}\) die Facette verdeckt (Ray-Cast). **Parallelisierbar pro Facette**.

**Regime-Adapter** (100–150 km):  
\[
\mathbf{F}=(1-\beta)\,\mathbf{F}_{\mathrm{FM}}+\beta\,\mathbf{F}_{\mathrm{corr}}(\mathrm{Kn});\quad
\beta=\exp(-\gamma\,\mathrm{Kn}),\;\gamma\sim\mathcal{O}(1)
\]
\(\mathbf{F}_{\mathrm{corr}}\) aus R13-Tabellen/Surrogats.

---

## 3) Software-Architektur (C++17/20)

**Module**
- `atm/` Atmosphären-Wrapper: NRLMSIS/DTM/JB2008 + HWM14; einheitliches Interface `AtmosphereState`.
- `gsi/` GSI-Kerne (Tabellen + Interpolator / kleines ONNX-NN).
- `geom/` Mesh-Import (OBJ/STL), Facettennormale, BVH; optional **Embree** (oder `nanort`) für Okklusion.
- `solver/` Panel-Integrator (seriell & parallel), Regime-Adapter, UQ (MC-Samples).
- `core/` Mathe/Typen (Eigen/handgerollt), Units (SI).
- `cli/` Kommandozeile (JSON-Input/Output), Benchmarks, Validierung.

**Datenflüsse**
```
Input: orbit_state, time, geom.mesh, materials.json, config.json
   -> AtmosphereState(h, lat, lon, t, F10.7, Kp)  --> {n_s, T, rho, wind}
   -> For each facet: view(θ_i), Ma, τ -> GSI Kern -> dF_i,s
   -> Occlusion (ray-cast) -> sum F, M
   -> Regime-Adapter (Kn) -> F_corr
   -> UQ (optional Monte Carlo)
Output: F, M, partials, variances
```

---

## 4) Parallelisierung & Performance

- **Task-Parallel** pro Facette (und/oder pro UQ-Sample).  
  - CPU: **OpenMP** `#pragma omp parallel for` (Chunking über Facetten).  
  - **Work-stealing** (TBB) optional.  
- **Vektorisiert** (Eigen/hand-SIMD) für Kern-Eval & Vektoralgebra.  
- **Okklusion**: BVH + **Embree** (single ray per facet) → sehr günstig; alternativ `nanort` (header-only).  
- **GPU (optional)**: CUDA-Kernel pro Facette (SoA), Textur-Lookup für Kerne; Host-Fallback beibehalten.

---

## 5) Validierung & Referenzen

**Vergleichs-Fälle**
- Flachplatte/Würfel im FM-Regime → Referenzkoeffizienten nach **Sentman**.  
- Sensitivität auf \(\alpha_n,\alpha_t,T_w\) vs. **CLL**-Vorhersagen.  
- Gesamtkräfte/Momente vs. publizierte Panel-Tools/NASA-Berichte (z. B. Aerodynamic Torques).  
- Dichte/Wind-Realismus mit **NRLMSIS 2.0 / DTM2020 / JB2008 / HWM14**.  
- DSMC-„Truth“ (nur als Vergleich): **Bird**-Standardreferenz.

---

## 6) Konfiguration (Beispiel)

```json
{
  "geometry": "satellite.stl",
  "cg": [0.0, 0.0, 0.0],
  "materials": {
    "default": {"alpha_n": 0.9, "alpha_t": 0.8, "Tw_K": 300.0}
  },
  "atmosphere": {
    "model": "NRLMSIS2",
    "winds": "HWM14",
    "indices": {"F10_7": 120.0, "Kp": 3}
  },
  "state": {"alt_km": 400.0, "lat_deg": 0.0, "lon_deg": 0.0, "utc": "2025-09-12T12:00:00Z",
            "V_sat_mps": [7500.0, 0.0, 0.0]},
  "solver": {"occlusion": "embree", "uq_samples": 0, "threads": "auto"}
}
```

---

## 7) Schritt-für-Schritt-Plan für den Code-Assistenten (C++)

> **Hinweis:** Bitte jeden Schritt vollständig umsetzen & testen, bevor der nächste folgt.

### Schritt 0 — Projekt-Skeleton
- Lege Repo an: `fmx/`  
- Unterordner: `core/ atm/ gsi/ geom/ solver/ cli/ tests/ cmake/ thirdparty/`  
- **CMake** (C++20), Abhängigkeiten optional: Eigen, OpenMP, Embree (optional), nlohmann/json, fmt.

### Schritt 1 — Kern-Typen & Utils
- `core/types.hpp`: `Vec3`, `Mat3`, `Aabb`, `Facet{A, n, r_center, material_id}`.  
- `core/units.hpp`: Konstanten (k_B, m_O, …) SI-konform.  
- `core/json.hpp`: Settings laden.

### Schritt 2 — Atmosphären-Interface
- `atm/Atmosphere.hpp` (reines Interface):  
  `AtmosphereState evaluate(alt_km, lat, lon, utc, F10_7, Kp)` → `{rho, T, wind, species[]}`.  
- `atm/NRLMSIS2.cpp` (Stub + klare TODO-Schnittstelle für spätere Kopplung).  
- `atm/HWM14.cpp` (Stub für Winde).

### Schritt 3 — Geometrie & Okklusion
- `geom/Mesh.cpp` (OBJ/STL-Loader; Recompute-Normals; Flächen; BVH).  
- `geom/Occluder.hpp`: Strategie-Interface (`none|bvh|embree`).  
- `geom/EmbreeOccluder.cpp` (optional; guard per CMake).

### Schritt 4 — GSI-Kerne
- `gsi/Sentman.hpp/.cpp`: direkte Formeln (Baseline).  
- `gsi/CLL.hpp/.cpp`: CLL-basierte Koeffizienten-Eval (Tabellen-Interpolation; Datenformat: HDF5).  
- `gsi/KernelSet.hpp`: Lookup: `(theta, Ma, tau, alpha_n, alpha_t, species)` → `(C_N,C_T, qdot)`.

### Schritt 5 — Panel-Solver (seriell)
- `solver/PanelSolver.hpp/.cpp`:  
  - pro Facette Sichtwinkel \(\theta\), Okklusions-Ray, Ma/τ berechnen, GSI-Kern abfragen, Summation \(\mathbf{F},\mathbf{M}\).  
  - Rückgaben inkl. Partials \(\partial \mathbf{F}/\partial \rho\).

### Schritt 6 — Regime-Adapter
- `solver/RegimeAdapter.hpp/.cpp`: Knudsen-Schätzung via mittlere freie Weglänge; Blend \( \beta(\mathrm{Kn}) \).  
- `gsi/R13Corr.hpp` (Stub): 1D-Tabellen für Korrekturen im Bereich Kn≈0.3–3.

### Schritt 7 — **Parallelisierung**
- Ersetze serielle Facetten-Schleifen durch **OpenMP** (oder TBB), reduziere auf \(\mathbf{F},\mathbf{M}\).  
- Parallelisiere **UQ-Samples** (äußere Schleife).  
- (Optional) SIMD-Friendly Arrays of Structs → Struct of Arrays.

### Schritt 8 — CLI & I/O
- `cli/fmx.cpp`: JSON einlesen, Solver aufrufen, JSON ausgeben (Kräfte/Momente, ggf. Varianzen).  
- `--bench` (Timer, Threads); `--validate <case>` (Standardfälle).

### Schritt 9 — Tests & Validierung
- Unit-Tests (Math, GSI-Kerne, Occlusion).  
- **Referenzfälle**: Flachplatte/Würfel (Sentman) & Publikationen (Torques).

### Schritt 10 — (Optional) GPU
- CUDA-Kernel pro Facette (SoA), Textur-Lookup für Kerne; Host-Fallback beibehalten.

---

## 8) Theoriedetails für die Implementierung

### 8.1 Sentman-Koeffizienten (Kurzform)
Für ein driftendes Maxwell-Gas liefern die Sentman-Integrale geschlossene Ausdrücke für Normal-/Scheranteile am Impulsfluss basierend auf \(\theta,\mathrm{Ma},\tau,\alpha\). Implementiere **numerisch stabile** Varianten (Vermeidung katastrophaler Auslöschung bei \(\theta\to\pi/2\)), prüfe Grenzfälle (specular/diffus).

### 8.2 CLL-Kern
Der **Cercignani–Lampis**-Kern modelliert die Verteilung reflektierter Geschwindigkeiten mit separaten normal/tangentialen Energie-Akk.-Koeffizienten; Implementierung via tabellierter Moment-Integrale/Werte → schnelle Interpolation.

### 8.3 Winde & Indizes
**HWM14** liefert \(\mathbf{u}_\text{wind}(h,\varphi,\lambda,t)\). Solar/Geomagnetik-Treiber: **F10.7** & **Kp/Ap**. Stelle die CLI so ein, dass Indizes direkt übergeben werden können (später: Online-Fetcher).

### 8.4 Unsicherheit
- **Dichte-Unsicherheit** als \(\delta\log\rho\sim\mathcal{N}(0,\sigma^2)\), \(\sigma\) modellabhängig (DTM2020-UQ als Referenzrahmen).  
- **Material/GSI**: Spannweiten für \(\alpha_n,\alpha_t\); **UQ** über MC-Samples.

---

## 9) Benchmarks (Zielwerte)

- **10⁴ Facetten**, 5 Spezies: **<10 ms** auf 8–16C CPU; **<1 ms** auf GPU (ohne Okklusion).  
- Okklusions-Ray pro Facette mit Embree: **~1–3 ms** (einmaliger BVH-Build amortisiert).

---

## 10) Beispiel-API (C++)

```cpp
struct Species { double rho, m; }; // kg/m3, kg
struct AtmosphereState {
  double T; Vec3 wind;
  std::array<Species, MAX_SPECIES> sp;
};

class Atmosphere {
public:
  virtual AtmosphereState evaluate(const GeoTime& gt, const Indices& idx) const = 0;
};

struct Material { double alpha_n, alpha_t, Tw; };

struct Facet { Vec3 n, r; double A; int mat_id; bool shadowed; };

struct Result { Vec3 F, M; Covariance cov; };

class PanelSolver {
public:
  Result solve(const Mesh& mesh, const AtmosphereState& atm,
               const Kinematics& kin, const Regime& reg, const UQ& uq) const;
};
```

---

## 11) Validierungs-Checkliste

- **Einfachkörper**: Flachplatte/Cube — stimmen \(C_D, C_L, C_M\) mit Sentman-/CLL-Referenzen?  
- **Okklusionstest**: Rotierender Würfel — Vergleich „mit/ohne Schatten“.  
- **Atmosphärenwechsel**: 200/400/600 km, Tag/Nacht, aktiv/ruhig (F10.7/Kp-Szenarien).  
- **Übergang**: 120–140 km, Blend-Funktion monotone Annäherung, keine Artefakte.  
- **Cross-Model**: NRLMSIS 2.0 vs. DTM2020 vs. JB2008.

---

## 12) Lizenz & Drittbibliotheken

- **Embree** (Apache-2.0), optional; **nanort** (MIT), optional.  
- Eigen (MPL2), nlohmann/json (MIT), fmt (MIT).  
- Diese Bibliothek: empfohlene Lizenz **BSD-3-Clause** oder **Apache-2.0**.

---

## 13) Literatur & Quellen (Auswahl)

**Atmosphäre & Winde**  
Emmert, J. T., et al. (2021). NRLMSIS 2.0: A Whole-Atmosphere Empirical Model. *JGR Space Physics*.  
Drob, D. P., et al. (2015). HWM14 – An update to the Horizontal Wind Model. *Earth and Space Science*.  
Bruinsma, S. (2015). DTM-2013 thermosphere model. *JSWSC*.  
Bruinsma, S., & Boniface, C. (2021). DTM2020 and uncertainties. *JSWSC*.  
Bowman, B. R., et al. (2008). JB2008 Thermospheric Density Model. *AIAA*.  

**GSI / Frei-molekular**  
Sentman, L. H. (1961). *Free-Molecule Flow Theory…*. LMSC-448514 / NASA.  
Cercignani, C., Lampis, M. (1972). Scattering kernels for gas–surface interactions. *TTSP*.  
Lord, R. (1991/95). Extensions to the CLL model. *Phys. Fluids A*.  
NASA CR-313 (1965). *Aerodynamic Torques on Satellites*.  

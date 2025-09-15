# FMX – Theorie & Methodik (100–600 km, frei‑molekular bis Übergang)

Diese Datei erklärt **Schritt für Schritt** die Theorie hinter der geplanten schnellen, parallelisierten Methode zur Berechnung atmosphärischer Kräfte/Momente auf Satelliten und begründet **warum** jeder Schritt nötig ist und **was er bringt**. Am Ende stehen vollständige Literaturangaben.

---

## 1) Umgebung: Thermosphäre & Treiber

**Was machen wir?**  
Wir bestimmen Temperatur, Dichte und Zusammensetzung der neutralen Atmosphäre (N₂, O₂, O, He, H …) sowie großskalige **Winde** entlang der Bahn (Zeit/Ort). Dafür nutzen wir etablierte **Empirikmodelle**: NRLMSIS 2.0 (T, Dichten), DTM 2020 (T, ρ, UQ), JB2008 (ρ) und **HWM14** (Winde). Solar-/geomagnetische Treiber sind **F10.7** (10.7 cm Radioflux) und **Kp/Ap**.  

**Warum?**  
Ohne realistische Randbedingungen sind Kräfte/Momente systematisch falsch. Die Modelle kondensieren Jahrzehnte an Messungen zu robusten Prädiktoren für 100–600 km.

**Was bringt es?**  
– Konsistente Inputs (T, ρ, nₛ, Winde) für alle Facetten/Spezies.  
– Einfache Variation nach Raumwetter (F10.7/Kp) und damit schnelle Sensitivitäts-/UQ‑Analysen.

**Kernformeln/Begriffe**  
Knudsen-Zahl \(\mathrm{Kn}=\lambda/L\) (frei‑molekular ab \(\mathrm{Kn}\gtrsim 10\)); mittlere freie Weglänge \(\lambda\) klassisch \(\lambda\propto (n\sigma)^{-1}\). Zusammensetzung: unter ~200 km N₂/O₂, darüber zunehmend O, >400 km He/H.  

---

## 2) Kinematik: Relativgeschwindigkeit & Spezies‑Machzahl

**Was machen wir?**  
Wir bilden die **Relativgeschwindigkeit** zwischen Satellit und Atmosphäre
\(\mathbf{c}=\mathbf{V}_\text{sat}-\mathbf{u}_\text{wind}\). Für jede Spezies \(s\) berechnen wir \(\mathrm{Ma}_s=\|\mathbf{c}\|/\sqrt{kT/m_s}\) und das Stagnationsimpulsmaß \(p_{\infty,s}=\rho_s\|\mathbf{c}\|^2\).

**Warum?**  
Kräfte im frei‑molekularen Regime skalieren mit \(\rho \|\mathbf{c}\|^2\) und hängen vom Winkel zwischen \(-\hat{\mathbf{c}}\) und der Flächennormale ab.

**Was bringt es?**  
– Richtige Spezies‑Gewichtung (O vs. N₂ etc.).  
– Korrekte Scher-/Normalanteile in den GSI‑Modellen.

---

## 3) Gas–Oberflächen‑Wechselwirkung (GSI)

**Was machen wir?**  
Für jede Facette mit Fläche \(A_i\), Normale \(\mathbf{n}_i\) und Wandtemperatur \(T_w\) integrieren wir den Molekülimpulsfluss mit **Sentman** (Basis) und optional **CLL** (Cercignani–Lampis–Lord) für realistischere Energie-/Impuls‑Akk. Die GSI liefert **Koeffizienten** \(C_N, C_T\) (normal/tangential) als Funktion von \(\theta\), \(\mathrm{Ma}_s\), \(\tau=T_w/T\) sowie \(\alpha_n,\alpha_t\) (Akk.-Koeffizienten).

\[
\mathrm{d}\mathbf{F}_{i,s}=p_{\infty,s}A_i\!\left[
C_{N,s}\,\mathbf{n}_i + C_{T,s}\,\hat{\mathbf{t}}_i\right],\qquad
\mathbf{F}=\sum_{i,s}\mathrm{d}\mathbf{F}_{i,s}
\]

**Warum?**  
Die GSI bestimmt Widerstand, Seitenkraft, Auftrieb und Wärmefluss. Sentman ist analytisch schnell; CLL bildet reale Oberflächen/Teil‑Diffuse/Teil‑Spekulare Effekte besser ab.

**Was bringt es?**  
– **Schnell** (Sentman geschl. Formeln) und **präzise** (CLL‑Tabellen/Surrogat).  
– Erweiterbar auf materialspezifische \(\alpha_n,\alpha_t\), \(T_w\) (z. B. AO‑erodierte Oberflächen).

---

## 4) Okklusion (Schatten) durch Geometrie

**Was machen wir?**  
Wir testen je Facette, ob der Strahl in Richtung \(-\hat{\mathbf{c}}\) die Facette trifft oder durch andere Geometrie abgeschattet wird (Ray‑Casting mit BVH/Embree). Abgeschattete Facetten tragen **nicht** bei.

**Warum?**  
Im frei‑molekularen Anströmen ist line‑of‑sight maßgeblich: verdeckte Flächen werden nicht von Molekülen getroffen.

**Was bringt es?**  
– Realistische Kräfte/Momente bei komplexen, konkaven Geometrien.  
– Sehr schneller Zusatzaufwand durch beschleunigte Strahl‑Geometrie‑Schnitt‑Tests.

---

## 5) Übergangsregime (≈100–150 km): R13‑Korrektur‑Blend

**Was machen wir?**  
Für \(\mathrm{Kn}\sim 0.1\ldots 10\) blenden wir zwischen frei‑molekularer Lösung und einer Korrektur aus **regularisierten 13‑Momenten (R13)**:
\(\mathbf{F}=(1-\beta)\mathbf{F}_\mathrm{FM}+\beta\,\mathbf{F}_\mathrm{corr}(\mathrm{Kn})\), \(\beta=\exp(-\gamma\,\mathrm{Kn})\).

**Warum?**  
Reine FM‑Theorie versagt bei dichterer Luft; R13 ist stabiler als Burnett und bildet Kn‑Effekte (Wärmefluss/Spannungen) konsistent ab.

**Was bringt es?**  
– Bessere Genauigkeit 100–150 km ohne teures DSMC.  
– Konsistenter Übergang ohne Artefakte.

---

## 6) Unsicherheiten (UQ)

**Was machen wir?**  
Wir ziehen **Ensembles** der Atmosphärenzustände (z. B. DTM2020‑UQ; MSIS‑UQ‑Kalibrierung) sowie Spannweiten für \(\alpha_n,\alpha_t,T_w\) und propagieren sie über Monte‑Carlo (parallel).

**Warum?**  
Thermosphäre ist hochvariabel; Dichte‑Bias/Scatter dominieren den Fehlerhaushalt vieler Missionen.

**Was bringt es?**  
– **Fehlerbalken** für Kräfte/Momente.  
– Robuste Design‑/Kontroll‑Entscheidungen (Worst‑Case, Prozentile).

---

## 7) Numerik & Parallelisierung

**Was machen wir?**  
– **Panelweise Parallelität** (OpenMP/TBB/CUDA).  
– **Lookup‑Kerne** für GSI (Tabellen/Surrogate), numerisch stabile Spezialfunktionen (erf, expm1) für große Ma und Grenzwinkel.  
– **Ray‑Casting** via Embree/nanort (BVH).

**Warum?**  
Leistung skaliert ~linear mit Facettenzahl; Ray‑Tests sind O(log N) pro Strahl; GSI‑Lookup vermeidet teure Integrale zur Laufzeit.

**Was bringt es?**  
– Millisekunden‑Antwortzeiten für 10⁴–10⁵ Facetten auf Mehrkern‑CPUs.  
– Optionale GPU‑Beschleunigung ohne Algorithmuswechsel.

---

## 8) Validierung

**Was machen wir?**  
– Flachplatte/Würfel (Sentman‑Referenz).  
– Literaturfälle (Aerodynamic Torques) & Panel‑Tools (z. B. ADBSat) vs. DSMC.  
– Cross‑Vergleich NRLMSIS 2.0 vs. DTM 2020 vs. JB2008; HWM14‑Winde.

**Warum?**  
Nachweis von Bias/Varianz und Gültigkeitsbereich (z. B. Probleme bei tiefen Konkavitäten).

**Was bringt es?**  
– Rückführbare Genauigkeit, reproduzierbare Benchmarks, dokumentierte Grenzen.

---

## 9) Wichtige Gleichungen (kompakt)

- \(\mathrm{Kn}=\lambda/L\), \(\lambda \approx \frac{kT}{\sqrt{2}\pi d^2 p}\).  
- \(\mathbf{c}=\mathbf{V}_\text{sat}-\mathbf{u}_\text{wind}\), \(\mathrm{Ma}_s=\|\mathbf{c}\|/\sqrt{kT/m_s}\).  
- Panelkraft je Spezies: \(\mathrm{d}\mathbf{F}_{i,s}= \rho_s\|\mathbf{c}\|^2A_i\,[C_{N,s}\mathbf{n}_i+C_{T,s}\hat{\mathbf{t}}_i]\).  
- Übergangs‑Blend: \(\mathbf{F}=(1-\beta)\mathbf{F}_\mathrm{FM}+\beta\,\mathbf{F}_\mathrm{corr}(\mathrm{Kn})\).

---

## Literatur

**Atmosphäre & Winde**  
Emmert, J.T., et al. (2021). *NRLMSIS 2.0: A Whole‑Atmosphere Empirical Model of Temperature and Neutral Species Densities*. JGR Space Physics.  
Bruinsma, S., Boniface, C. (2021). *The operational and research DTM‑2020 thermosphere models*. JSWSC.  
Bowman, B.R., et al. (2008). *JB2008 Thermospheric Density Model*. AIAA 2008‑6438 / SET TR‑2008‑002.  
Drob, D.P., et al. (2015). *An update to the Horizontal Wind Model: HWM14*. Earth and Space Science.  
Tapping, K.F. (2013). *The 10.7 cm solar radio flux (F10.7)*. Space Weather.  
NOAA SWPC. *The K‑index* (Tech Note); NCEI Space‑Weather Produkte.

**GSI / Frei‑molekular**  
Sentman, L.H. (1961). *Free Molecule Flow Theory and Its Application to the Determination of Aerodynamic Forces*. LMSC‑448514.  
Cercignani, C., Lampis, M. (1972). *Scattering kernels for gas–surface interactions*. Transport Theory & Statistical Physics.  
Lord, R.G. (1991, 1995). *Extensions of the Cercignani–Lampis kernel*. Physics of Fluids A.  
NASA CR‑313 (1965). *Aerodynamic Torques on Satellites*.  
Sharipov, F. (2002, 2003). *Applications of the CLL boundary condition*. Eur. J. Mech. B/Fluids.

**Übergang / R13 / DSMC**  
Struchtrup, H.; Torrilhon, M. (2003, 2004). *Regularization of Grad’s 13‑moment equations*; *Shock structure with R13*. Phys. Fluids; J. Fluid Mech.  
Struchtrup, H. (2005). *Macroscopic Transport Equations for Rarefied Gas Flows*. Springer.  
Bird, G.A. (1994/1998). *Molecular Gas Dynamics and the Direct Simulation of Gas Flows*. OUP.

**Materialeffekte (AO)**  
Banks, B.A., et al. (2003/2004). *Atomic Oxygen Effects/Interactions in LEO*. NASA TM.  
de Groh, K.K., et al. (2019). *Atomic Oxygen Erosion Data From the MISSE 2–8 Missions*. NASA.

**Okklusion / Ray‑Casting**  
Thomas, R.E. (1993). *Mutual Molecular Shadowing with Sentman’s Equations*. NASA.  
Intel Embree. *High‑Performance Ray Tracing Kernels* (Apache‑2.0).  
nanort (MIT). *Header‑only BVH kernel*.

**Panel‑Methoden‑Validierung**  
Sinpetru, L.A., et al. (2021, 2022). *ADBSat: Methodology & Validation*. arXiv/Acta Astronautica.

---

## Hinweise zur Implementierung

- Numerische Stabilität: erf/expm1‑Varianten, Grenzfälle \(\theta\to 0,\pi/2\).  
- SoA‑Datenlayout für SIMD/GPU; Reduktionen mit Kahan/compensated summation.  
- UQ: Log‑Normal für Dichte‑Scatter; Percentile‑Reporting (P5/50/95).

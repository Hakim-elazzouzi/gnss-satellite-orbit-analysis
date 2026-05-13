# 🛰️ Project 2 — GNSS Satellite Orbit Analysis

> **Keplerian Elements · Ground Tracks · Orbital Speed · RAAN Wheel · Constellation Geometry | 2026-01-01**

---

## 📌 Overview

Building on Project 1 (satellite XYZ positions), this project performs a full
**orbital mechanics analysis** of the GPS constellation. Instead of just computing
where the satellites are, we now analyse *how* they move and *why* the constellation
is designed the way it is.

---

## 📐 Key Orbital Mechanics

### The 6 Keplerian elements:

| Element | Symbol | Typical GPS value |
|---------|--------|-------------------|
| Semi-major axis | a | ~26,560 km |
| Eccentricity | e | ~0.01 (nearly circular) |
| Inclination | i | ~55° |
| RAAN | Ω₀ | 0°, 60°, 120°, 180°, 240°, 300° (6 planes) |
| Argument of perigee | ω | varies |
| Mean anomaly | M₀ | varies |

### Key formulas:
```
Orbital speed   v = √(μ/a)           ≈ 3.874 km/s
Orbital period  T = 2π√(a³/μ)        ≈ 11h 58m
Ground track    λ = atan2(Y, X)       longitude
                φ = atan2(Z, √X²+Y²)  latitude
```

---

## 🖼️ Output Plots

### Plot 1 — Keplerian Elements Dashboard
Four bar charts (a, e, i, Ω₀) for all 32 GPS satellites.
Shows how uniformly the constellation maintains its design parameters.

### Plot 2 — Ground Tracks
Sub-satellite points for all GPS satellites over 24 hours on a world map.
Sinusoidal paths — bounded by ±55° latitude (the orbital inclination).
Auckland station marked as a reference point.

### Plot 3 — Orbital Speed and Period
Bar charts of v and T per satellite. All should be nearly identical (~3.874 km/s,
~11.97 h) — any outlier indicates a satellite in a different slot.

### Plot 4 — RAAN Wheel (Polar Plot)
Polar plot showing each satellite's RAAN. The 6-plane Walker structure
(A through F, 60° apart) is immediately visible.

### Plot 5 — Constellation Geometry
- Scatter: Inclination vs RAAN (coloured by semi-major axis)
- Histogram: Eccentricity distribution across the full fleet

---

## 📂 File Structure

```
gnss-satellite-orbit-analysis/
│
├── output/
│   ├── plot1_keplerian_elements.png
│   ├── plot2_ground_tracks.png
│   ├── plot3_speed_period.png
│   ├── plot4_raan_wheel.png
│   └── plot5_constellation_geometry.png
├── src/
│   └── project2_gnss_satellites_orbit_analysis.py     ← Main python
├── requirements.txt
├── LICENSE
└── README.md
```

---

## ⚙️ How to Run

### 1. Install dependencies
```bash
pip install -r requirements.txt
```

### 2. Set your navigation file path
In **Step 2** of the notebook:
```python
nav_path = "../../data/brdc0010.26n"
```

### 3. Run all cells
```bash
jupyter notebook src/project2_orbit_analysis.ipynb
```

---

## 🛠️ Dependencies

| Package | Purpose |
|---------|---------|
| `georinex` | Parse RINEX 2/3 navigation files |
| `numpy` | Numerical computations |
| `pandas` | Time series management |
| `matplotlib` | 2D and polar publication-quality plotting |

---

## 👤 Author

**Hakim El Azzouzi**  
MSc Global Navigation Satellite Systems  
Mohammed First University, Oujda, Morocco  
📧 elazzouzihakim10@gmail.com  
🔗 [linkedin.com/in/Hakim-El-Azzouzi](https://linkedin.com/in/Hakim-El-Azzouzi)  
📍 Luxembourg 🇱🇺

---

## 📜 License

MIT License — see [LICENSE](LICENSE) for details.

---

## 🔗 GNSS Navigation RINEX Series

| # | Project |
|---|---------|
| 1 | GPS Satellite Position Computation |
| **2** | **GNSS Satellite Orbit Analysis** ← You are here |
| 3 | GNSS Azimuth & Elevation Analysis |
| 4 | GNSS Navigation Data Pre-processing |
| 5 | GNSS Single-Point Positioning |

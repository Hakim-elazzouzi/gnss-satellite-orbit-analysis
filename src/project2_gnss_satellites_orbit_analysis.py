=============================================================================
🛰️ Project 2 — GNSS Satellite Orbit Analysis
=============================================================================
 Author   : Hakim El Azzouzi
 Degree   : MSc Global Navigation Satellite Systems
            Mohammed First University, Oujda, Morocco
 Email    : elazzouzihakim10@gmail.com
 LinkedIn : https://linkedin.com/in/Hakim-El-Azzouzi
 Location : Luxembourg 🇱🇺
-----------------------------------------------------------------------------
 File type    : RINEX 2 GPS Navigation
 Date         : 2026-01-01 (Day of Year 001)
 Satellites   : GPS only (G01–G32)
 Records      : ~400 ephemeris sets (each valid ~2 hours)
 Source       : IGS CDDIS broadcast combined file
-----------------------------------------------------------------------------
 Description
 -----------
In Project 1 we computed the raw ECEF X, Y, Z positions of every GPS satellite.
This project goes deeper: we analyse the **orbital mechanics** of those trajectories
and produce a set of diagnostic plots.

**5 plots produced:**

| Plot | What It Shows |
|------|---------------|
| 🌍 Orbital plane geometry | Inclination, RAAN, and orbital plane orientation for each satellite |
| 📐 Keplerian elements | How key orbital parameters vary across the constellation |
| 🔄 Ground track | Sub-satellite point latitude vs longitude on a world map |
| ⏱️ Orbital period & velocity | Speed and period computed from the ephemeris |
| 📊 Constellation coverage | How many satellites are above 10° elevation from Auckland |

---

## 📐 Orbital Mechanics Background

A GPS satellite orbit is described by **6 Keplerian elements**:

| Element | Symbol | Meaning | Typical GPS value |
|---------|--------|---------|-------------------|
| Semi-major axis | a | Size of the orbit | ~26,560 km |
| Eccentricity | e | How elliptical | ~0.01 (nearly circular) |
| Inclination | i | Tilt vs equatorial plane | ~55° |
| RAAN | Ω₀ | Longitude of ascending node | varies (6 planes × 60°) |
| Arg. of perigee | ω | Orientation within orbital plane | varies |
| Mean anomaly | M₀ | Position at reference epoch | varies |

GPS uses **6 orbital planes** separated by 60° in RAAN, each containing
~4 satellites — giving at least 4 satellites visible from anywhere on Earth at all times.

### Orbital velocity and period:
```
v = √(μ/a)   ≈  3.874 km/s     orbital speed
T = 2π√(a³/μ) ≈  11h 58m        sidereal orbital period
```

### Ground track:
The sub-satellite point (SSP) is the point on Earth directly below the satellite.
Because Earth rotates underneath the orbit, the SSP traces a sinusoidal curve on a map.
GPS satellites repeat their ground track every **sidereal day** (≈ 23h 56m 4s).
---
-----------------------------------------------------------------------------
 **About the projects**
 ----------------------
# Step1: Install & Import Libraries
# Step2: Load the RINEX File
# Step3: Compute Quality Metrics for Every Satellite
# Step4: Plot 1: Coverage & SNR Dashboard
# Step5: Plot 2: Data Gap Map
# Step6: Plot 3: Quality Score Summary Scatter
# Step7: Generate the Text Quality Report
=============================================================================
# ───────────────────────────────────
# Step 1 — Install & Import Libraries
# ───────────────────────────────────
# Uncomment if running for the first time:
# !pip install --upgrade georinex numpy matplotlib pandas

import georinex as gr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as mticker
import warnings

warnings.filterwarnings('ignore')
plt.rcParams.update({'figure.dpi': 120})

# ── GPS Physical Constants (ICD-GPS-200) ───────────────────────────────────
MU      = 3.986005e14           # Earth gravitational constant [m³/s²]
OMEGA_E = 7.2921151467e-5       # Earth rotation rate [rad/s]
R_EARTH = 6_371_000.0           # Mean Earth radius [m]

# Derived GPS nominal values (for reference)
A_GPS_NOM = 26_560_000.0        # Nominal GPS semi-major axis [m]
V_GPS_NOM = np.sqrt(MU / A_GPS_NOM)           # Nominal orbital speed [m/s]
T_GPS_NOM = 2 * np.pi * np.sqrt(A_GPS_NOM**3 / MU)  # Nominal period [s]

print('✅ Libraries loaded')
print(f'   Nominal GPS orbital radius : {A_GPS_NOM/1e3:.0f} km')
print(f'   Nominal GPS orbital speed  : {V_GPS_NOM/1e3:.3f} km/s')
print(f'   Nominal GPS orbital period : {T_GPS_NOM/3600:.4f} hours  (~11h 58m)')

# ──────────────────────────────────────────────────────────────────
# Step 2 — Load the Navigation File and Reuse the Position Algorithm
# ──────────────────────────────────────────────────────────────────

# The navigation file path

nav_path = "/brdc0010.26n"   # ← change this if needed

print("⏳ Loading navigation file...")
nav = gr.load(nav_path)
print("✅ Done!")
print()

all_sv   = nav.sv.values
gps_sats = sorted([s for s in all_sv if s.startswith('G')])
print(f"GPS satellites found: {len(gps_sats)}")
print("  " + "  ".join(gps_sats))


# ── ICD-GPS-200 position algorithm (same as Project 1) ─────────────────────
def gps_satellite_position(eph, t_gps):
    sqrtA    = float(eph['sqrtA']);       a        = sqrtA ** 2
    e        = float(eph['Eccentricity'])
    M0       = float(eph['M0'])
    delta_n  = float(eph['DeltaN'])
    omega    = float(eph['omega'])
    i0       = float(eph['Io'])
    IDOT     = float(eph['IDOT'])
    Omega0   = float(eph['Omega0'])
    OmegaDot = float(eph['OmegaDot'])
    toe      = float(eph['Toe'])
    Cuc = float(eph['Cuc']); Cus = float(eph['Cus'])
    Crc = float(eph['Crc']); Crs = float(eph['Crs'])
    Cic = float(eph['Cic']); Cis = float(eph['Cis'])

    tk = t_gps - toe
    if tk >  302400: tk -= 604800
    if tk < -302400: tk += 604800

    n  = np.sqrt(MU / a**3) + delta_n
    Mk = M0 + n * tk

    Ek = Mk
    for _ in range(50):
        Ek_new = Mk + e * np.sin(Ek)
        if abs(Ek_new - Ek) < 1e-12: break
        Ek = Ek_new
    Ek = Ek_new

    vk    = np.arctan2(np.sqrt(1-e**2)*np.sin(Ek), np.cos(Ek)-e)
    Phi_k = vk + omega

    uk = Phi_k + Cus*np.sin(2*Phi_k) + Cuc*np.cos(2*Phi_k)
    rk = a*(1-e*np.cos(Ek)) + Crs*np.sin(2*Phi_k) + Crc*np.cos(2*Phi_k)
    ik = i0 + IDOT*tk + Cis*np.sin(2*Phi_k) + Cic*np.cos(2*Phi_k)

    xk = rk * np.cos(uk)
    yk = rk * np.sin(uk)

    Ok = Omega0 + (OmegaDot - OMEGA_E)*tk - OMEGA_E*toe

    X = xk*np.cos(Ok) - yk*np.cos(ik)*np.sin(Ok)
    Y = xk*np.sin(Ok) + yk*np.cos(ik)*np.cos(Ok)
    Z = yk*np.sin(ik)
    return X, Y, Z


def find_best_ephemeris(nav, sat, t_gps):
    try:
        sat_data = nav.sel(sv=sat)
    except Exception:
        return None
    toes     = sat_data['Toe'].values
    best_idx = np.nanargmin(np.abs(toes - t_gps))
    if np.abs(toes[best_idx] - t_gps) > 7200:
        return None
    return sat_data.isel(time=best_idx)


print()
print("✅ Position algorithm ready")

# ─────────────────────────────────────────────────
# Step 3 — Extract Keplerian Elements per Satellite
# ─────────────────────────────────────────────────
# For each satellite, read the mean Keplerian elements from all ephemeris records
keplerian = {}   # sat → dict of orbital parameters

print(f"{'Sat':<6} {'a [km]':>10} {'e':>8} {'i [°]':>8} {'Ω₀ [°]':>10} {'ω [°]':>10} {'T [h]':>8} {'v [km/s]':>10}")
print("-" * 75)

for sat in gps_sats:
    try:
        sd = nav.sel(sv=sat)
    except Exception:
        continue

    # Semi-major axis from sqrtA
    sqrtA_vals = sd['sqrtA'].values
    sqrtA_vals = sqrtA_vals[~np.isnan(sqrtA_vals)]
    if len(sqrtA_vals) == 0:
        continue
    a = float(np.mean(sqrtA_vals))**2   # [m]

    # Eccentricity
    e_vals = sd['Eccentricity'].values
    e = float(np.nanmean(e_vals))

    # Inclination
    i_vals = sd['Io'].values
    inc    = float(np.nanmean(i_vals))   # [rad]

    # Right ascension of ascending node
    O_vals = sd['Omega0'].values
    Omega0 = float(np.nanmean(O_vals))   # [rad]

    # Argument of perigee
    w_vals = sd['omega'].values
    omega  = float(np.nanmean(w_vals))   # [rad]

    # Derived quantities
    T = 2 * np.pi * np.sqrt(a**3 / MU)   # orbital period [s]
    v = np.sqrt(MU / a)                   # mean orbital speed [m/s]

    keplerian[sat] = {
        'a':      a,
        'e':      e,
        'inc':    inc,
        'Omega0': Omega0,
        'omega':  omega,
        'T':      T,
        'v':      v,
    }

    print(f"  {sat:<4}  {a/1e3:>10.2f}  {e:>8.5f}  {np.degrees(inc):>8.3f}"
          f"  {np.degrees(Omega0) % 360:>10.3f}  {np.degrees(omega) % 360:>10.3f}"
          f"  {T/3600:>8.4f}  {v/1e3:>10.4f}")

print()
print(f"✅ Keplerian elements extracted for {len(keplerian)} satellites")

# ────────────────────────────────────────────────────────────
# Step 4 — Plot 1: Keplerian Elements Across the Constellation
# ────────────────────────────────────────────────────────────
sat_names = list(keplerian.keys())
n = len(sat_names)
x = np.arange(n)

palette = plt.cm.tab20(np.linspace(0, 1, n))

fig, axes = plt.subplots(2, 2, figsize=(18, 10), facecolor='#0d1117')
fig.suptitle(
    'GPS Keplerian Orbital Elements — Full Constellation\n'
    'Extracted from Broadcast Ephemeris | 2026-01-01',
    fontsize=14, fontweight='bold', color='#ffffff'
)

panels = [
    ('a',      'Semi-Major Axis [km]',       1e3,  26400, 26700,
     'All satellites orbit at nearly the same altitude (~26,560 km)'),
    ('e',      'Eccentricity [-]',            1.0,  0.0,   0.025,
     'Very small eccentricity — GPS orbits are nearly perfectly circular'),
    ('inc',    'Inclination [°]',            np.pi/180, 50, 60,
     'All planes inclined at ~55° to the equatorial plane'),
    ('Omega0', 'RAAN Ω₀ [°]',               np.pi/180, 0,  365,
     'RAAN spans 0°–360°, spread across 6 planes separated by 60°'),
]

for ax, (key, ylabel, scale, ymin, ymax, note) in zip(axes.flat, panels):
    ax.set_facecolor('#111827')
    ax.tick_params(colors='#aaaaaa', labelsize=8)
    ax.grid(True, axis='y', color='#222222', linewidth=0.5)
    for spine in ax.spines.values():
        spine.set_edgecolor('#333333')

    vals = [keplerian[s][key] / scale for s in sat_names]

    # Special handling: wrap RAAN to 0–360°
    if key == 'Omega0':
        vals = [v % 360 for v in vals]

    bars = ax.bar(x, vals,
                  color=[palette[i] for i in range(n)],
                  edgecolor='#0d1117', linewidth=0.6)

    ax.set_xticks(x)
    ax.set_xticklabels(sat_names, rotation=60, ha='right', fontsize=7)
    for tick, color in zip(ax.get_xticklabels(), palette):
        tick.set_color(color)

    ax.set_ylabel(ylabel, color='#aaaaaa', fontsize=10)
    ax.set_ylim(ymin, ymax)
    ax.set_title(note, color='white', fontsize=9, pad=4)

    # Reference line for semi-major axis and inclination
    if key == 'a':
        ax.axhline(26560, color='#FFEB3B', ls='--', lw=1.2,
                   label='Nominal 26,560 km')
        legend = ax.legend(fontsize=8, framealpha=0.3,
                           facecolor='#1a1a2e', edgecolor='#444444')
        for t in legend.get_texts(): t.set_color('white')
    if key == 'inc':
        ax.axhline(55.0, color='#FFEB3B', ls='--', lw=1.2,
                   label='Design inclination 55°')
        legend = ax.legend(fontsize=8, framealpha=0.3,
                           facecolor='#1a1a2e', edgecolor='#444444')
        for t in legend.get_texts(): t.set_color('white')

plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig('../content/output/plot1_keplerian_elements.png',
            dpi=150, bbox_inches='tight', facecolor=fig.get_facecolor())
plt.show()

print('✅ Saved: ../content/output/plot1_keplerian_elements.png')
print()
print('💡 Interpretation:')
print('   • Semi-major axis: all ~26,560 km — GPS design requirement')
print('   • Eccentricity: all < 0.02 — nearly circular, by design')
print('   • Inclination: all ~55° — allows coverage up to 55° latitude')
print('   • RAAN: spread every ~30–60° across 6 planes — ensures global coverage')

# ──────────────────────────────────
# Step 5 — Plot 2: GPS Ground Tracks
# ──────────────────────────────────
"""The sub-satellite point (SSP) is the point
on Earth's surface directly below each satellite.
As the satellite orbits and Earth rotates underneath it,
the SSP traces a sinusoidal path on a world map.

We convert ECEF (X, Y, Z) to geodetic latitude and longitude:

longitude λ = atan2(Y, X)           [degrees]
latitude  φ = atan2(Z, √(X²+Y²))   [degrees, approximate spherical]

"""
# Time grid: every 15 minutes over 24 hours
GPS_WEEK_START = 345600.0
STEP_SEC       = 900.0
N_EPOCHS       = 96

t_gps_epochs   = GPS_WEEK_START + np.arange(N_EPOCHS) * STEP_SEC
utc_timestamps = pd.date_range('2026-01-01 00:00', periods=N_EPOCHS, freq='15min')

# Compute XYZ and convert to lat/lon
print("⏳ Computing ground tracks...")

ground_tracks = {}   # sat → {'lat': array, 'lon': array}

for sat in gps_sats:
    lats = np.full(N_EPOCHS, np.nan)
    lons = np.full(N_EPOCHS, np.nan)

    for i, t in enumerate(t_gps_epochs):
        eph = find_best_ephemeris(nav, sat, t)
        if eph is None:
            continue
        try:
            X, Y, Z = gps_satellite_position(eph, t)
            # Convert ECEF → spherical geodetic (approximate, ignoring Earth flattening)
            lons[i] = np.degrees(np.arctan2(Y, X))          # longitude [°]
            lats[i] = np.degrees(np.arctan2(Z, np.sqrt(X**2 + Y**2)))  # latitude [°]
        except Exception:
            pass

    ground_tracks[sat] = {'lat': lats, 'lon': lons}

print(f"✅ Ground tracks computed for {len(ground_tracks)} satellites")

# ── Plot ────────────────────────────────────────────────────────────────────
palette_gt = plt.cm.plasma(np.linspace(0.1, 0.9, len(gps_sats)))

fig, ax = plt.subplots(figsize=(18, 9), facecolor='#0d1117')
ax.set_facecolor('#0d1117')

# World map background grid
ax.set_xlim(-180, 180)
ax.set_ylim(-90, 90)
ax.axhline(0, color='#1a4a6b', lw=0.6, ls='--', alpha=0.6)   # equator
ax.axhline(+55, color='#2d6a4f', lw=0.5, ls=':', alpha=0.5)   # GPS inclination
ax.axhline(-55, color='#2d6a4f', lw=0.5, ls=':', alpha=0.5)

for lon_line in range(-180, 181, 30):
    ax.axvline(lon_line, color='#1a1a2e', lw=0.5, alpha=0.6)
for lat_line in range(-90, 91, 30):
    ax.axhline(lat_line, color='#1a1a2e', lw=0.5, alpha=0.6)

# Plot each satellite ground track
for i, sat in enumerate(gps_sats):
    lats = ground_tracks[sat]['lat']
    lons = ground_tracks[sat]['lon']
    mask = ~np.isnan(lats)

    # Plot the track — break where there are large longitude jumps (antimeridian crossing)
    lon_diff = np.abs(np.diff(lons[mask]))
    breaks   = np.where(lon_diff > 180)[0] + 1

    segments_lon = np.split(lons[mask], breaks)
    segments_lat = np.split(lats[mask], breaks)

    for seg_lon, seg_lat in zip(segments_lon, segments_lat):
        ax.plot(seg_lon, seg_lat,
                color=palette_gt[i], lw=1.0, alpha=0.75)

    # Mark starting position (00:00 UTC)
    if not np.isnan(lats[0]):
        ax.scatter(lons[0], lats[0],
                   color=palette_gt[i], s=25, zorder=5, edgecolors='white', lw=0.3)
        ax.text(lons[0] + 1, lats[0] + 1, sat[1:],
                fontsize=5.5, color=palette_gt[i], alpha=0.9)

# Mark Auckland station
ax.scatter(174.76, -36.88, color='#FFEB3B', s=80, marker='*',
           zorder=10, label='AUCK00NZL (Auckland)')
ax.text(176, -35.5, 'Auckland', color='#FFEB3B', fontsize=9, fontweight='bold')

# Labels
ax.set_xlabel('Longitude [°]', color='#aaaaaa', fontsize=11)
ax.set_ylabel('Latitude [°]', color='#aaaaaa', fontsize=11)
ax.tick_params(colors='#aaaaaa')
ax.set_xticks(range(-180, 181, 30))
ax.set_yticks(range(-90, 91, 30))
for spine in ax.spines.values():
    spine.set_edgecolor('#333333')

ax.set_title(
    'GPS Satellite Ground Tracks — Sub-Satellite Points\n'
    '2026-01-01 | 15-minute intervals | Dots = position at 00:00 UTC\n'
    'Green dotted lines = ±55° inclination limit',
    fontsize=12, fontweight='bold', color='#ffffff'
)

legend = ax.legend(fontsize=10, loc='lower left',
                   framealpha=0.4, facecolor='#1a1a2e', edgecolor='#444444')
for t in legend.get_texts(): t.set_color('white')

plt.tight_layout()
plt.savefig('../content/output/plot2_ground_tracks.png',
            dpi=150, bbox_inches='tight', facecolor=fig.get_facecolor())
plt.show()

print('✅ Saved: output/plot2_ground_tracks.png')
print()
print('💡 Interpretation:')
print('   • Sinusoidal shape = orbit inclined at 55° + Earth rotation underneath')
print('   • Latitude never exceeds ±55° — GPS cannot fly over the poles')
print('   • Satellites repeat their ground track every sidereal day (~23h 56m)')
print('   • Multiple satellites visible from Auckland (★) at all times')

# ───────────────────────────────────────────────────────
# Step 6 — Plot 3: Orbital Speed and Period per Satellite
# ───────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(16, 6), facecolor='#0d1117')
fig.suptitle(
    'GPS Orbital Speed and Period per Satellite\n'
    'Computed from Broadcast Ephemeris | 2026-01-01',
    fontsize=13, fontweight='bold', color='#ffffff'
)

sat_names_k = list(keplerian.keys())
x_k = np.arange(len(sat_names_k))
palette_k = plt.cm.tab20(np.linspace(0, 1, len(sat_names_k)))

for ax in axes:
    ax.set_facecolor('#111827')
    ax.tick_params(colors='#aaaaaa', labelsize=8)
    ax.grid(True, axis='y', color='#222222', linewidth=0.5)
    for spine in ax.spines.values():
        spine.set_edgecolor('#333333')

# Left: orbital speed
speeds = [keplerian[s]['v'] / 1e3 for s in sat_names_k]   # km/s
axes[0].bar(x_k, speeds,
            color=[palette_k[i] for i in range(len(sat_names_k))],
            edgecolor='#0d1117', linewidth=0.6)
axes[0].axhline(V_GPS_NOM / 1e3, color='#FFEB3B', ls='--', lw=1.5,
                label=f'Nominal v = {V_GPS_NOM/1e3:.3f} km/s')
axes[0].set_xticks(x_k)
axes[0].set_xticklabels(sat_names_k, rotation=60, ha='right', fontsize=7)
for tick, color in zip(axes[0].get_xticklabels(), palette_k):
    tick.set_color(color)
axes[0].set_ylabel('Orbital Speed [km/s]', color='#aaaaaa', fontsize=11)
axes[0].set_title('v = √(μ/a)  —  depends only on semi-major axis',
                  color='white', fontsize=10)
axes[0].set_ylim(3.8, 3.95)

legend0 = axes[0].legend(fontsize=9, framealpha=0.3,
                          facecolor='#1a1a2e', edgecolor='#444444')
for t in legend0.get_texts(): t.set_color('white')

# Right: orbital period
periods = [keplerian[s]['T'] / 3600 for s in sat_names_k]   # hours
axes[1].bar(x_k, periods,
            color=[palette_k[i] for i in range(len(sat_names_k))],
            edgecolor='#0d1117', linewidth=0.6)
axes[1].axhline(T_GPS_NOM / 3600, color='#FFEB3B', ls='--', lw=1.5,
                label=f'Nominal T = {T_GPS_NOM/3600:.4f} h')
axes[1].set_xticks(x_k)
axes[1].set_xticklabels(sat_names_k, rotation=60, ha='right', fontsize=7)
for tick, color in zip(axes[1].get_xticklabels(), palette_k):
    tick.set_color(color)
axes[1].set_ylabel('Orbital Period [hours]', color='#aaaaaa', fontsize=11)
axes[1].set_title('T = 2π√(a³/μ)  —  ~11h 58m  (exactly half a sidereal day)',
                  color='white', fontsize=10)
axes[1].set_ylim(11.9, 12.1)

legend1 = axes[1].legend(fontsize=9, framealpha=0.3,
                          facecolor='#1a1a2e', edgecolor='#444444')
for t in legend1.get_texts(): t.set_color('white')

plt.tight_layout(rect=[0, 0, 1, 0.94])
plt.savefig('../content/output/plot3_speed_period.png',
            dpi=150, bbox_inches='tight', facecolor=fig.get_facecolor())
plt.show()

print('✅ Saved: output/plot3_speed_period.png')
print()
print('📊 Speed and period summary:')
print(f"   Nominal speed  : {V_GPS_NOM/1e3:.4f} km/s")
print(f"   Actual range   : {min(speeds):.4f} – {max(speeds):.4f} km/s")
print(f"   Nominal period : {T_GPS_NOM/3600:.4f} h")
print(f"   Actual range   : {min(periods):.4f} – {max(periods):.4f} h")
print()
print('💡 The tiny spread (~0.001 km/s) shows all satellites are maintained')
print('   in very tightly controlled orbital slots by ground control.')

# ────────────────────────────────────────────────────────
# Step 7 — Plot 4: Orbital Plane Distribution (RAAN Wheel)
# ────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 10),
                       subplot_kw={'projection': 'polar'},
                       facecolor='#0d1117')
ax.set_facecolor('#0d1117')

palette_w = plt.cm.tab20(np.linspace(0, 1, len(sat_names_k)))

# GPS plane colours (A–F)
PLANE_COLORS = {
    'A': '#2196F3', 'B': '#F44336', 'C': '#4CAF50',
    'D': '#FF9800', 'E': '#9C27B0', 'F': '#00BCD4',
}

for i, sat in enumerate(sat_names_k):
    raan_deg = np.degrees(keplerian[sat]['Omega0']) % 360
    raan_rad = np.radians(raan_deg)
    r        = keplerian[sat]['a'] / 1e6   # Mm, used as radial distance for display

    ax.scatter(raan_rad, r,
               color=palette_w[i], s=120, zorder=5,
               edgecolors='white', linewidths=0.5)
    ax.text(raan_rad, r + 0.12, sat,
            ha='center', va='center',
            fontsize=7, color=palette_w[i])

# Reference circles at nominal GPS radius
r_nom = A_GPS_NOM / 1e6
theta_circle = np.linspace(0, 2*np.pi, 300)
ax.plot(theta_circle, np.full(300, r_nom),
        color='#FFEB3B', lw=0.8, ls='--', alpha=0.5)

# Mark the 6 nominal GPS planes every 60°
for plane_idx, letter in enumerate('ABCDEF'):
    angle = np.radians(plane_idx * 60)
    ax.plot([angle, angle], [0, r_nom + 0.5],
            color='#333333', lw=0.8, ls=':')
    ax.text(angle, r_nom + 0.7, f'Plane {letter}',
            ha='center', va='center',
            color='#888888', fontsize=8)

# Styling
ax.set_theta_zero_location('N')     # 0° at top (north)
ax.set_theta_direction(-1)          # clockwise (like a compass)
ax.set_rlabel_position(90)
ax.tick_params(colors='#aaaaaa', labelsize=8)
ax.set_rticks([24, 25, 26, 26.56, 27, 28])
ax.set_rlim(0, 29)
ax.grid(color='#222222', linewidth=0.5)

# Custom angle labels
ax.set_thetagrids(range(0, 360, 30),
                  [f'{a}°' for a in range(0, 360, 30)],
                  color='#888888', fontsize=8)

ax.set_title(
    'GPS Orbital Plane Distribution — RAAN Wheel\n'
    'Radial axis = semi-major axis [Mm]  |  Angular axis = RAAN Ω₀\n'
    'Six nominal planes (A–F) separated by 60°',
    fontsize=12, fontweight='bold', color='#ffffff', pad=20
)

plt.savefig('../content/output/plot4_raan_wheel.png',
            dpi=150, bbox_inches='tight', facecolor=fig.get_facecolor())
plt.show()

print('✅ Saved: output/plot4_raan_wheel.png')
print()
print('💡 Interpretation:')
print('   • Satellites cluster in 6 groups ~60° apart in RAAN')
print('   • This is the GPS Walker constellation design')
print('   • Each cluster = one orbital plane (A through F)')
print('   • Gaps would indicate a missing satellite or repositioning manoeuvre')

# ──────────────────────────────────────────────────
# Step 8 — Plot 5: Constellation Geometry Statistics
# ──────────────────────────────────────────────────

fig, axes = plt.subplots(1, 2, figsize=(16, 6), facecolor='#0d1117')
fig.suptitle(
    'GPS Constellation Geometry Statistics\n'
    'Orbital plane orientation parameters | 2026-01-01',
    fontsize=13, fontweight='bold', color='#ffffff'
)

for ax in axes:
    ax.set_facecolor('#111827')
    ax.tick_params(colors='#aaaaaa')
    ax.grid(True, color='#222222', linewidth=0.5)
    for spine in ax.spines.values():
        spine.set_edgecolor('#333333')

# Left: Inclination vs RAAN scatter
raan_vals = [np.degrees(keplerian[s]['Omega0']) % 360 for s in sat_names_k]
inc_vals  = [np.degrees(keplerian[s]['inc'])           for s in sat_names_k]
a_vals    = [keplerian[s]['a'] / 1e3                   for s in sat_names_k]

sc = axes[0].scatter(raan_vals, inc_vals,
                     c=a_vals, cmap='plasma',
                     s=120, edgecolors='white', linewidths=0.5, zorder=5)

for sat, raan, inc in zip(sat_names_k, raan_vals, inc_vals):
    axes[0].annotate(sat, (raan, inc),
                     xytext=(4, 4), textcoords='offset points',
                     fontsize=6.5, color='#aaaaaa')

cbar = plt.colorbar(sc, ax=axes[0], pad=0.02)
cbar.set_label('Semi-major axis [km]', color='#e0e0e0')
plt.setp(cbar.ax.get_yticklabels(), color='#aaaaaa')
cbar.ax.yaxis.set_tick_params(color='#aaaaaa')

axes[0].axhline(55.0, color='#FFEB3B', ls='--', lw=1.0,
                alpha=0.7, label='Design inclination 55°')
axes[0].set_xlabel('RAAN Ω₀ [°]', color='#aaaaaa', fontsize=11)
axes[0].set_ylabel('Inclination i [°]', color='#aaaaaa', fontsize=11)
axes[0].set_title('Inclination vs RAAN\n'
                  'Colour = semi-major axis — all should be ~26,560 km',
                  color='white', fontsize=10)
axes[0].set_xlim(-5, 365)
axes[0].set_ylim(50, 60)

legend0 = axes[0].legend(fontsize=9, framealpha=0.3,
                          facecolor='#1a1a2e', edgecolor='#444444')
for t in legend0.get_texts(): t.set_color('white')

# Right: Eccentricity distribution
ecc_vals = [keplerian[s]['e'] for s in sat_names_k]
bins     = np.linspace(0, 0.03, 20)
axes[1].hist(ecc_vals, bins=bins,
             color='#2196F3', edgecolor='#0d1117', alpha=0.85)
axes[1].axvline(np.mean(ecc_vals), color='#FFEB3B', ls='--', lw=1.5,
                label=f'Mean e = {np.mean(ecc_vals):.5f}')
axes[1].set_xlabel('Eccentricity e [-]', color='#aaaaaa', fontsize=11)
axes[1].set_ylabel('Number of satellites', color='#aaaaaa', fontsize=11)
axes[1].set_title('Eccentricity distribution\n'
                  'All values near zero — GPS orbits are nearly circular',
                  color='white', fontsize=10)

legend1 = axes[1].legend(fontsize=9, framealpha=0.3,
                          facecolor='#1a1a2e', edgecolor='#444444')
for t in legend1.get_texts(): t.set_color('white')

plt.tight_layout(rect=[0, 0, 1, 0.94])
plt.savefig('../content/output/plot5_constellation_geometry.png',
            dpi=150, bbox_inches='tight', facecolor=fig.get_facecolor())
plt.show()

print('✅ Saved: output/plot5_constellation_geometry.png')
print()
print('📊 Constellation geometry summary:')
print(f"   Inclination  : {min(inc_vals):.3f}° – {max(inc_vals):.3f}°  (mean {np.mean(inc_vals):.3f}°)")
print(f"   RAAN         : {min(raan_vals):.1f}° – {max(raan_vals):.1f}°")
print(f"   Eccentricity : {min(ecc_vals):.6f} – {max(ecc_vals):.6f}  (mean {np.mean(ecc_vals):.6f})")
print(f"   Semi-major a : {min(a_vals):.2f} – {max(a_vals):.2f} km")


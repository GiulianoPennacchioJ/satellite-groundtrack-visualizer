from datetime import datetime, timedelta
import math
from pathlib import Path

# Earth parameters (WGS-84)
MU_EARTH_KM3_S2 = 398600.4418      # km^3/s^2
EARTH_RADIUS_KM = 6378.137         # km

# I/O directories
BASE_DIR = Path(__file__).resolve().parent
IN_DIR = BASE_DIR / "IN"
OUT_DIR = BASE_DIR / "OUT"
# Also look for shared TLEs in the sibling ConstellationVisibility project.
TLE_SEARCH_DIRS = [
    IN_DIR,
    BASE_DIR / ".." / "ConstellationVisibility-local" / "IN",
]
OUT_DIR.mkdir(parents=True, exist_ok=True)


def analysis_anchor_datetime():
    """
    Reproducible anchor: fixed hardcoded date at 00:00 UTC.
    """
    return datetime(2025, 11, 26, 0, 0, 0)


# ------------------------------------------------------------
# 1) Read TLEs (3-line format: NAME, LINE1, LINE2)
# ------------------------------------------------------------
def read_tles(filename):
    """
    Read a text file with TLEs in 3-line format:
        NAME
        LINE1
        LINE2
    Returns a list of dicts: [{ 'name','line1','line2' }, ...]
    """
    sats = []
    with open(filename, "r") as f:
        lines = [l.strip() for l in f if l.strip()]

    i = 0
    while i < len(lines):
        line = lines[i]
        # Allow also 2-line TLEs (no name)
        if line.startswith("1 ") or line.startswith("2 "):
            name = f"SAT_{len(sats)+1}"
            line1 = lines[i]
            line2 = lines[i + 1]
            i += 2
        else:
            name = line
            line1 = lines[i + 1]
            line2 = lines[i + 2]
            i += 3

        sats.append({"name": name, "line1": line1, "line2": line2})

    return sats


def resolve_tle_path(filename):
    """
    Find the first existing TLE file across the known input directories.
    """
    filename = Path(filename)
    if filename.is_absolute():
        return filename if filename.exists() else None

    for base in TLE_SEARCH_DIRS:
        candidate = base / filename
        if candidate.exists():
            return candidate
    return None


# ------------------------------------------------------------
# 2) Parse TLE → Keplerian elements at epoch
# ------------------------------------------------------------
def parse_tle_epoch(line1):
    """
    Parse TLE epoch (line 1) to a datetime in UTC.
    """
    year_str = line1[18:20]
    doy_str = line1[20:32]

    year = int(year_str)
    # TLE convention: 57-99 → 19xx, 00-56 → 20xx
    if year < 57:
        year += 2000
    else:
        year += 1900

    doy = float(doy_str)
    day_int = int(doy)
    frac_day = doy - day_int

    epoch = datetime(year, 1, 1) + timedelta(days=day_int - 1 + frac_day)
    return epoch


def parse_tle_elements(tle):
    """
    From raw TLE lines, extract:
      a, e, i, RAAN, arg_perigee, M0, n, epoch
    Angles in radians, a in km, n in rad/s.
    """
    line1 = tle["line1"]
    line2 = tle["line2"]

    epoch = parse_tle_epoch(line1)

    # Slices follow standard TLE column layout
    i_deg      = float(line2[8:16])
    raan_deg   = float(line2[17:25])
    e_str      = line2[26:33]
    argp_deg   = float(line2[34:42])
    M_deg      = float(line2[43:51])
    n_rev_day  = float(line2[52:63])

    # Eccentricity has implied decimal point
    e = float("0." + e_str.strip())

    i_rad     = math.radians(i_deg)
    raan_rad  = math.radians(raan_deg)
    argp_rad  = math.radians(argp_deg)
    M0_rad    = math.radians(M_deg)

    # Mean motion [rev/day] → [rad/s]
    n_rad_s = 2.0 * math.pi * n_rev_day / 86400.0

    # Kepler's 3rd law: n^2 a^3 = mu
    a = (MU_EARTH_KM3_S2 / (n_rad_s**2)) ** (1.0 / 3.0)

    return {
        "name": tle["name"],
        "epoch": epoch,
        "a_km": a,
        "e": e,
        "i_rad": i_rad,
        "raan_rad": raan_rad,
        "argp_rad": argp_rad,
        "M0_rad": M0_rad,
        "n_rad_s": n_rad_s,
    }


# ------------------------------------------------------------
# 3) Simple ECI → ECEF → geodetic helpers
# ------------------------------------------------------------
def gmst_angle(dt):
    """
    Approximate Greenwich Mean Sidereal Time angle [rad].
    Good enough for visualization.
    """
    # Julian Date (simple formula)
    jd = (367 * dt.year
          - int((7 * (dt.year + int((dt.month + 9) / 12))) / 4)
          + int((275 * dt.month) / 9)
          + dt.day
          + 1721013.5
          + (dt.hour + dt.minute / 60.0 + dt.second / 3600.0) / 24.0)

    T = (jd - 2451545.0) / 36525.0

    gmst_sec = (67310.54841
                + (876600.0 * 3600 + 8640184.812866) * T
                + 0.093104 * T * T
                - 6.2e-6 * T * T * T)

    gmst_sec = gmst_sec % 86400.0
    gmst_rad = (gmst_sec / 240.0) * math.pi / 180.0
    return gmst_rad


def eci_to_ecef(r_eci_km, dt):
    """
    Rotate position from ECI to ECEF using GMST.
    r_eci_km: [x, y, z] in km
    """
    theta = gmst_angle(dt)
    cos_t = math.cos(theta)
    sin_t = math.sin(theta)

    x, y, z = r_eci_km

    x_ecef = cos_t * x + sin_t * y
    y_ecef = -sin_t * x + cos_t * y
    z_ecef = z

    return [x_ecef, y_ecef, z_ecef]


def ecef_to_geodetic(r_ecef_km):
    """
    Very simple spherical conversion: ECEF [km] → (lat, lon, alt_km).
    Good enough for a workshop visualization.
    """
    x, y, z = r_ecef_km

    r_xy = math.hypot(x, y)
    r = math.sqrt(x * x + y * y + z * z)

    if r_xy == 0.0:
        lat = math.copysign(math.pi / 2, z)
        lon = 0.0
    else:
        lat = math.atan2(z, r_xy)
        lon = math.atan2(y, x)

    alt_km = r - EARTH_RADIUS_KM

    lat_deg = math.degrees(lat)
    lon_deg = math.degrees(lon)
    return lat_deg, lon_deg, alt_km


# ------------------------------------------------------------
# 4) Kepler propagation
# ------------------------------------------------------------
def solve_kepler(M, e, tol=1e-8, max_iter=20):
    """
    Solve Kepler's equation M = E - e sin(E) for E.
    Newton-Raphson method.
    """
    E = M
    for _ in range(max_iter):
        f = E - e * math.sin(E) - M
        fp = 1 - e * math.cos(E)
        if fp == 0:
            break
        dE = f / fp
        E_new = E - dE
        if abs(dE) < tol:
            return E_new
        E = E_new
    return E


def propagate_satellite(tle, start_dt, end_dt, step_seconds=60):
    """
    Propagate one satellite (2-body Kepler) from start_dt to end_dt.
    Returns list of dicts with lat/lon/alt over time.
    """
    elems = parse_tle_elements(tle)

    epoch   = elems["epoch"]
    a       = elems["a_km"]
    e       = elems["e"]
    i_rad   = elems["i_rad"]
    raan    = elems["raan_rad"]
    argp    = elems["argp_rad"]
    M0      = elems["M0_rad"]
    n_rad_s = elems["n_rad_s"]

    points = []
    t = start_dt

    while t <= end_dt:
        # 1) Mean anomaly at time t
        dt_sec = (t - epoch).total_seconds()
        M = M0 + n_rad_s * dt_sec
        M = math.fmod(M, 2.0 * math.pi)
        if M < 0:
            M += 2.0 * math.pi

        # 2) Solve Kepler
        E = solve_kepler(M, e)
        cosE = math.cos(E)

        # 3) Radius and true anomaly
        r_km = a * (1 - e * cosE)
        cos_v = (cosE - e) / (1 - e * cosE)
        sin_v = (math.sqrt(1 - e * e) * math.sin(E)) / (1 - e * cosE)
        v = math.atan2(sin_v, cos_v)

        # 4) Position in ECI
        u = argp + v
        cos_u, sin_u = math.cos(u), math.sin(u)
        cos_i, sin_i = math.cos(i_rad), math.sin(i_rad)
        cos_O, sin_O = math.cos(raan), math.sin(raan)

        x = r_km * (cos_O * cos_u - sin_O * sin_u * cos_i)
        y = r_km * (sin_O * cos_u + cos_O * sin_u * cos_i)
        z = r_km * (sin_u * sin_i)

        r_eci = [x, y, z]

        # Speed magnitude from vis-viva (km/s)
        speed_km_s = math.sqrt(MU_EARTH_KM3_S2 * (2.0 / r_km - 1.0 / a))

        # Inertial (ECI) spherical lat/lon (no Earth rotation)
        eci_lat_deg, eci_lon_deg, eci_alt_km = ecef_to_geodetic(r_eci)

        # 5) ECI → ECEF → geodetic
        r_ecef = eci_to_ecef(r_eci, t)
        lat_deg, lon_deg, alt_km = ecef_to_geodetic(r_ecef)

        points.append({
            "name": elems["name"],
            "time": t,
            "eci_x_km": x,
            "eci_y_km": y,
            "eci_z_km": z,
            "eci_lat_deg": eci_lat_deg,
            "eci_lon_deg": eci_lon_deg,
            "eci_alt_km": eci_alt_km,
            "speed_km_s": speed_km_s,
            "lat_deg": lat_deg,
            "lon_deg": lon_deg,
            "alt_km": alt_km,
        })

        t += timedelta(seconds=step_seconds)

    return points


# ------------------------------------------------------------
# 5) Coastline loader and ground track plot
# ------------------------------------------------------------
def load_world_map_segments(filename=None):
    """
    Read coastline segments from a text file with 'lat lon' per line and blank
    lines separating segments. Returns [(lons, lats), ...].
    """
    if filename is None:
        filename = IN_DIR / "world_map.txt"
    segments = []
    lons, lats = [], []

    try:
        with open(filename, "r") as f:
            for line in f:
                stripped = line.strip()
                if not stripped:
                    if lons:
                        segments.append((lons, lats))
                        lons, lats = [], []
                    continue
                parts = stripped.split()
                if len(parts) != 2:
                    continue
                lat, lon = float(parts[0]), float(parts[1])
                lons.append(lon)
                lats.append(lat)
    except FileNotFoundError:
        return []

    if lons:
        segments.append((lons, lats))
    return segments


def plot_ground_tracks(points_by_sat, world_map_file=None, png_out=None):
    """
    Plot lon/lat ground tracks for each satellite using coastline data from
    world_map.txt. If png_out is provided, saves the figure instead of showing it.
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is required for plotting ground tracks.")
        return

    coastline_segments = load_world_map_segments(world_map_file)
    if not coastline_segments:
        print(f"Could not load coastline data from {world_map_file}.")

    fig, ax = plt.subplots(figsize=(11, 6), constrained_layout=True)
    ax.grid(linewidth=1.0)
    ax.set_axisbelow(True)

    for lon_seg, lat_seg in coastline_segments:
        ax.plot(lon_seg, lat_seg, color="#909090", linewidth=0.6,
                alpha=0.8, zorder=1)

    handles = []
    cmap = plt.get_cmap("tab10")
    sat_names = list(points_by_sat.keys())
    denom = max(1, len(sat_names) - 1)

    for idx, sat_name in enumerate(sat_names):
        pts = points_by_sat[sat_name]
        if not pts:
            continue
        lons = [p["lon_deg"] for p in pts]
        lats = [p["lat_deg"] for p in pts]
        color = cmap(idx / float(denom))
        h = ax.scatter(lons, lats, s=3, marker="o", linewidths=0,
                       alpha=0.85, zorder=2, label=sat_name, color=color)
        handles.append(h)

    ax.set_xlabel("Longitude [°]")
    ax.set_ylabel("Latitude [°]")
    ax.set_title("Ground Tracks")
    ax.set_xlim([-180, 180])
    ax.set_ylim([-90, 90])
    ax.legend()

    if png_out:
        fig.savefig(png_out, dpi=200, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show()


# ------------------------------------------------------------
# 6) Velocity time series plot
# ------------------------------------------------------------
def plot_velocity_timeseries(points_by_sat, png_out="velocity_timeseries.png"):
    """
    Plot speed [km/s] vs time for each satellite on separate subplots.
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is required for plotting velocity time series.")
        return

    sat_names = [k for k, v in points_by_sat.items() if v]
    if not sat_names:
        print("No points available for velocity plotting.")
        return

    fig, axes = plt.subplots(
        len(sat_names), 1, figsize=(12, 2.4 * len(sat_names)), sharex=True,
        constrained_layout=True
    )
    if len(sat_names) == 1:
        axes = [axes]

    for ax, sat_name in zip(axes, sat_names):
        pts = points_by_sat[sat_name]
        times = [p["time"] for p in pts]
        speeds = [p["speed_km_s"] for p in pts]
        ax.plot(times, speeds, color="#d1495b", linewidth=1.2)
        ax.set_ylabel("km/s")
        ax.set_title(sat_name, fontsize="small")
        ax.grid(True, linewidth=0.8, alpha=0.6)
        ax.ticklabel_format(style="plain", axis="y", useOffset=False)

    axes[-1].set_xlabel("UTC Time")
    fig.savefig(png_out, dpi=200, bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------------
# 7) Orbit report
# ------------------------------------------------------------
def orbit_summary(tle):
    """
    Compute orbital period and perigee/apogee altitudes (km) from TLE.
    """
    elems = parse_tle_elements(tle)
    a = elems["a_km"]
    e = elems["e"]
    n = elems["n_rad_s"]
    period_s = 2.0 * math.pi / n
    perigee_km = a * (1.0 - e) - EARTH_RADIUS_KM
    apogee_km = a * (1.0 + e) - EARTH_RADIUS_KM
    return {
        "name": elems["name"],
        "period_s": period_s,
        "perigee_km": perigee_km,
        "apogee_km": apogee_km,
        "a_km": a,
        "e": e,
        "i_deg": math.degrees(elems["i_rad"]),
        "raan_deg": math.degrees(elems["raan_rad"]),
        "argp_deg": math.degrees(elems["argp_rad"]),
        "M0_deg": math.degrees(elems["M0_rad"]),
    }


def write_orbit_report(summaries, filename=None):
    """
    Write orbital period (hours) and perigee/apogee altitudes for each satellite
    plus a second table with Keplerian elements.
    """
    if filename is None:
        filename = OUT_DIR / "orbit_report.txt"
    header = "Satellite Orbit Report"
    columns = "Name                     Period [h]   Perigee Alt [km]   Apogee Alt [km]"
    underline = "-" * len(columns)
    lines = [header, underline, columns, underline]

    for s in summaries:
        period_h = s["period_s"] / 3600.0
        lines.append(
            f"{s['name']:<24}{period_h:11.3f}{s['perigee_km']:18.2f}{s['apogee_km']:17.2f}"
        )

    # Second table: Keplerian elements at epoch
    lines.append("")
    kep_columns = ("Name                     a [km]        e        i [deg]    "
                   "RAAN [deg]  argPeri [deg]  M0 [deg]")
    kep_underline = "-" * len(kep_columns)
    lines.extend(["Keplerian Elements (at TLE epoch)", kep_underline, kep_columns, kep_underline])
    for s in summaries:
        lines.append(
            f"{s['name']:<24}"
            f"{s['a_km']:10.1f}"
            f"{s['e']:11.7f}"
            f"{s['i_deg']:12.3f}"
            f"{s['raan_deg']:12.3f}"
            f"{s['argp_deg']:14.3f}"
            f"{s['M0_deg']:10.3f}"
        )
    with open(filename, "w") as f:
        f.write("\n".join(lines) + "\n")


# ------------------------------------------------------------
# 8) KML writer
# ------------------------------------------------------------
def write_kml(points_by_sat, filename):
    """
    points_by_sat:
      {
        "SAT_NAME_1": [(lon_deg, lat_deg, alt_m), ...],
        "SAT_NAME_2": [...],
      }
    """
    header = """<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2">
<Document>
    <name>Satellite Orbits (Kepler)</name>
"""
    footer = """</Document>
</kml>
"""

    placemarks = []
    for sat_name, coords in points_by_sat.items():
        coord_strings = []
        for lon_deg, lat_deg, alt_m in coords:
            coord_strings.append(f"{lon_deg},{lat_deg},{alt_m}")
        coord_block = " ".join(coord_strings)

        placemark = f"""
    <Placemark>
        <name>{sat_name}</name>
        <Style>
            <LineStyle>
                <width>2</width>
            </LineStyle>
        </Style>
        <LineString>
            <tessellate>1</tessellate>
            <altitudeMode>absolute</altitudeMode>
            <coordinates>
                {coord_block}
            </coordinates>
        </LineString>
    </Placemark>
"""
        placemarks.append(placemark)

    with open(filename, "w") as f:
        f.write(header)
        for p in placemarks:
            f.write(p)
        f.write(footer)


# ------------------------------------------------------------
# 7) Example main
# ------------------------------------------------------------
if __name__ == "__main__":
    tle_file = resolve_tle_path("tle.txt")  # put your TLEs here (3-line format)

    if tle_file is None:
        search_paths = [str(path / "tle.txt") for path in TLE_SEARCH_DIRS]
        joined = "\n  - ".join(search_paths)
        print("Could not find tle.txt in any of:")
        print(f"  - {joined}")
        print("Put your TLEs there (3-line format).")
    else:
        try:
            sats = read_tles(tle_file)
        except FileNotFoundError:
            print(f"Could not find {tle_file}. Put your TLEs there (3-line format).")
        else:
            # Time window for propagation
            start = analysis_anchor_datetime()
            end = start + timedelta(hours=24)
            step_seconds = 5

            points_by_sat = {}
            summaries = []
            for tle in sats:
                pts = propagate_satellite(tle, start, end, step_seconds=step_seconds)
                points_by_sat[tle["name"]] = pts
                summaries.append(orbit_summary(tle))

            out_file = OUT_DIR / "orbits_kepler.kml"
            kml_tracks = {
                sat_name: [
                    (p["lon_deg"], p["lat_deg"], p["alt_km"] * 1000.0) for p in pts
                ]
                for sat_name, pts in points_by_sat.items()
            }

            write_kml(kml_tracks, out_file)
            print(f"Wrote {out_file}. Open it in Google Earth.")

            eci_out_file = OUT_DIR / "orbits_eci.kml"
            eci_kml_tracks = {
                sat_name: [
                    (p["eci_lon_deg"], p["eci_lat_deg"], p["eci_alt_km"] * 1000.0)
                    for p in pts
                ]
                for sat_name, pts in points_by_sat.items()
            }
            write_kml(eci_kml_tracks, eci_out_file)
            print(f"Wrote {eci_out_file} (tracks in ECI frame).")

            write_orbit_report(summaries, filename=OUT_DIR / "orbit_report.txt")
            print("Wrote OUT/orbit_report.txt.")

            plot_ground_tracks(
                points_by_sat,
                world_map_file=IN_DIR / "world_map.txt",
                png_out=OUT_DIR / "ground_tracks.png",
            )

            plot_velocity_timeseries(
                points_by_sat,
                png_out=OUT_DIR / "velocity_timeseries.png",
            )

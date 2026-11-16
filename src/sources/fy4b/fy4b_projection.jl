"""
    latlon_to_scan_angles(lat, lon, sat_lon=DEFAULT_SAT_LON; height=0.0)

Convert latitude/longitude to scan angles (x, y in radians) for geostationary projection.
Uses standard geostationary projection formula (similar to GOES-R ABI).
Returns: (scan_x, scan_y) in radians, or (NaN, NaN) if not visible.

`height` (metres) lifts the point along the local vertical before projecting, which gives where
a cloud top at that height over the spot appears in the image: the parallax-displaced position.
Over Hubei a 10 km top is displaced about 10 km toward the NW from 133E and about 8 km toward
the N from 105E. `height == 0` is the surface point, computed exactly as before.
"""
function latlon_to_scan_angles(lat, lon, sat_lon=DEFAULT_SAT_LON; height=0.0)
    # Convert to radians
    φ = deg2rad(lat)
    λ = deg2rad(lon - sat_lon)  # Relative longitude

    # WGS84 parameters
    re = EARTH_RADIUS           # Equatorial radius
    rp = 6356.75231414e3        # Polar radius (m)
    f = (re - rp) / re          # Flattening
    e2 = f * (2 - f)            # Eccentricity squared

    # Geocentric latitude (accounting for Earth flattening)
    φ_c = atan((1 - e2) * tan(φ))

    # Distance from Earth center to surface point
    rc = rp / sqrt(1 - e2 * cos(φ_c)^2)

    # Cartesian coordinates of surface point (Earth-centered, x-axis through 0°E)
    xs = rc * cos(φ_c) * cos(λ)
    ys = rc * cos(φ_c) * sin(λ)
    zs = rc * sin(φ_c)

    # Lift along the geodetic normal; skipped at zero so the surface path stays bit-identical.
    if height != 0
        xs += height * cos(φ) * cos(λ)
        ys += height * cos(φ) * sin(λ)
        zs += height * sin(φ)
    end

    # Visibility check: point must be on Earth side facing the satellite
    # The satellite is at (SAT_HEIGHT, 0, 0) in this coordinate system
    # Point is visible if it's not behind the Earth's limb
    dot_product = xs * (xs - SAT_HEIGHT) + ys^2 + zs^2
    if dot_product > 0
        # Check if behind Earth
        # Actually, we need to check if the line from satellite to point intersects Earth
        # Simplified: check the angle from satellite
        d2 = (SAT_HEIGHT - xs)^2 + ys^2 + zs^2
        cos_angle = (SAT_HEIGHT - xs) / sqrt(d2)
        if cos_angle < 0.156  # cos(81°) - maximum view angle
            return (NaN, NaN)
        end
    end

    # Scan angles from satellite perspective
    # x (E-W, west negative): angle in the equatorial plane
    scan_x = atan(ys / (SAT_HEIGHT - xs))

    # y (N-S, north positive): angle perpendicular to equatorial plane
    scan_y = atan(zs / sqrt((SAT_HEIGHT - xs)^2 + ys^2))

    return (scan_x, scan_y)
end

"""
    scan_angles_to_grid(scan_x, scan_y)

Convert scan angles (radians) to FY4B grid indices (0-based).
The grid is centered at the satellite sub-point (nadir).
"""
function scan_angles_to_grid(scan_x, scan_y)
    if isnan(scan_x) || isnan(scan_y)
        return (NaN, NaN)
    end

    # Convert scan angle to linear distance at nadir
    # At nadir: angle = distance / (H - R)
    # scan_angle = distance / SAT_DISTANCE
    # distance = scan_angle * SAT_DISTANCE

    # Grid spacing in radians (at nadir, RESOLUTION meters = scan_angle * SAT_DISTANCE)
    grid_spacing_rad = RESOLUTION / SAT_DISTANCE

    # Grid center
    center = (GRID_SIZE - 1) / 2.0

    # Convert scan angles to grid indices
    # Negative scan_x (west of nadir) -> smaller x index
    # Positive scan_y (north of nadir) -> smaller y index (image coordinates)
    x_idx = center + scan_x / grid_spacing_rad
    y_idx = center - scan_y / grid_spacing_rad

    return (x_idx, y_idx)
end

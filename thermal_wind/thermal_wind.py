"""Thermal wind exercise code for Mathematical Modelling of Geophysical
Fluids MPE2013 Workshop at African Institute for Mathematical Sciences.
"""
import matplotlib.pyplot as plt
import numpy as np


def thermal_wind(
    stn1_profile,
    stn2_profile,
    no_motion_depth=None,
    surface_delta=None,
):
    """Return the profile of horizontal velocity across the line between 2
    density profile stations by integrating the density profiles and using
    the specified boundary condition (level of no motion, or surface height
    difference) to calculate the velocity.

    :arg stn1_profile: Name of file containing the density profile at station 1
    :type stn1_profile: string

    :arg stn2_profile: Name of file containing the density profile at station 2
    :type stn2_profile: string

    :arg no_motion_depth: Depth at which level of no motion boundary condition
                          occurs
    :type no_motion_depth: number

    :arg surface_delta: Surface height difference between the 2 station
                        boundary condition
    :type surface_delta: number

    :returns: :obj:`(depth, v_vel)` Profile of horizontal velocity across
              the line between the 2 stations
    :rtype: tuple of :class:`numpy.ndarray`

    :raises: :exc:`IOError` if :obj:`stn1_profile` or :obj:`stn2_profile`
             file cannot be read

    :raises: :exc:`ValueError` if neither :obj:`no_motion_depth` or
             :obj:`surface_delta` are specified

    :raises: :exc:`ValueError` if both :obj:`no_motion_depth` and
             :obj:`surface_delta` are specified

    :raises: :exc:`ValueError` if depth intervals in the 2 station density
             profiles are unequal

    :raises: :exc:`ValueError` if :obj:`no_motion_depth` exceeds depth of
             shallowest density profile
    """
    # Geographical constants
    distance = 36.89e3  # distance between stations in m
    avg_latitude = 26  # average latitude of stations in degrees
    # Physical constants
    acc_grav = 9.8  # m / s^2
    rho0 = 1024  # kg / m^3
    Omega = 2 * np.pi / (60 * 60 * 24)
    fo = 2 * Omega * np.sin(np.pi * avg_latitude / 180)
    # Read the density profiles
    depth1, rho1 = read_density_profile(stn1_profile)
    depth2, rho2 = read_density_profile(stn2_profile)
    # Validate the argument values, read the density profiles, and calculate
    # depth to which they will be integrated
    v_vel_depth = validate_inputs(
        depth1, depth2, no_motion_depth, surface_delta)
    # Integrate the density from the surface to the depth of the shallowest
    # profile
    del_rho = rho2[:v_vel_depth.size] - rho1[:v_vel_depth.size]
    delta_density = np.empty_like(v_vel_depth)
    delta_density[0] = del_rho[0] * v_vel_depth[:1] / 2
    delta_density[1:-1] = (
        del_rho[1:-1] * (v_vel_depth[2:] - v_vel_depth[:-2]) / 2)
    delta_density[-1] = del_rho[-1] * (v_vel_depth[-1] - v_vel_depth[-2])
    sum_delta_density = delta_density.cumsum()
    # Apply the boundary condition
    if no_motion_depth is None:
        surface_effect = acc_grav * surface_delta / (fo * distance)
        no_motion_effect = 0
    else:
        mask = v_vel_depth >= no_motion_depth
        i_below = -v_vel_depth[mask].size
        i_above = i_below + 1
        proportion = (
            (no_motion_depth - v_vel_depth[i_above])
            / (v_vel_depth[i_below] - v_vel_depth[i_above]))
        sum_delta_density_no_motion = (
            sum_delta_density[i_above] + proportion
            * (sum_delta_density[i_below] - sum_delta_density[i_above]))
        no_motion_effect = (
            -acc_grav / (rho0 * fo) * sum_delta_density_no_motion / distance)
        surface_effect = 0
    # Calculate the horizontal velocity values
    v_vel = (
        surface_effect + no_motion_effect
        + acc_grav / (rho0 * fo) * sum_delta_density / distance)
    return v_vel_depth, v_vel


def validate_inputs(depth1, depth2, no_motion_depth, surface_delta):
    """Validate the inputs of the :func:`thermal_wind` function.

    * Exactly 1 of :obj:`no_motion_depth` or :obj:`surface_delta` must be
      specified

    * Depth intervals must be the same in both profiles

    * :obj`no_motion_depth` must not exceed the depth of the shallowest
      profile

    :arg depth1: Depths of the station 1 density profile
    :type depth1: :class:`numpy.ndarray`

    :arg depth2: Depths of the station 2 density profile
    :type depth2: :class:`numpy.ndarray`

    :arg no_motion_depth: Depth at which level of no motion boundary condition
                          occurs
    :type no_motion_depth: number

    :arg surface_delta: Surface height difference between the 2 station
                        boundary condition
    :type surface_delta: number

    :returns max_depth: Average of the deepest 2 levels in the shallowest
                        depth profile
    :rtype: float

    :raises: :exc:`ValueError` if neither :obj:`no_motion_depth` or
             :obj:`surface_delta` are specified

    :raises: :exc:`ValueError` if both :obj:`no_motion_depth` and
             :obj:`surface_delta` are specified

    :raises: :exc:`ValueError` if depth intervals in the 2 station density
             profiles are unequal

    :raises: :exc:`ValueError` if :obj:`no_motion_depth` exceeds depth of
             shallowest density profile
    """
    if (no_motion_depth is None and surface_delta is None
            or no_motion_depth is not None and surface_delta is not None):
        raise ValueError(
            'Specify either no motion depth or surface height difference')
    try:
        if depth1.size >= depth2.size:
            depth = depth2
            if not all(depth1[:depth2.size] == depth2):
                raise ValueError
        else:
            depth = depth1
            if not all(depth2[:depth1.size] == depth1):
                raise ValueError
    except ValueError:
        raise ValueError('Depth intervals must be the same in both profiles')
    if no_motion_depth > depth[-1]:
        raise ValueError('No motion depth is deeper than data')
    return depth


def read_density_profile(filename):
    """Return depth and density arrays read from filename.
    """
    depth = []
    density = []
    with open(filename) as f:
        for line in interesting_lines(f):
            deep, rho = map(float, line.split())
            depth.append(deep)
            density.append(rho)
    return np.array(depth), np.array(density)


def interesting_lines(f):
    for line in f:
        if line and not line.startswith('#'):
            yield line


def plot_velocity_profile(depth, velocity):
    """Plot the specified velocity component profile.

    :arg depth: Depths
    :type depth: :class:`numpy.ndarray`

    :arg velocity: Velocity component values
    :type velocity: :class:`numpy.ndarray`
    """
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(velocity, -depth)
    ax.set_xlabel('v [m/s]')
    ax.set_ylabel('z [m]')

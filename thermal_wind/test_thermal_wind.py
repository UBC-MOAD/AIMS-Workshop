"""Unit tests for thermal_wind module.
"""
import nose.tools
import numpy as np
import thermal_wind


@nose.tools.raises(ValueError)
def test_either_no_motion_depth_or_surface_delta():
    thermal_wind.validate_inputs(
        'foo', 'bar', no_motion_depth=280, surface_delta=0.2)


@nose.tools.raises(ValueError)
def test_one_of_no_motion_depth_or_surface_delta():
    depth1, rho1 = thermal_wind.read_density_profile('s109.dens')
    depth2, rho2 = thermal_wind.read_density_profile('s105.dens')
    thermal_wind.validate_inputs(
        depth1, depth2, no_motion_depth=None, surface_delta=None)


@nose.tools.raises(ValueError)
def test_unequal_depth_steps():
    depth1, rho1 = thermal_wind.read_density_profile('s109.dens')
    depth2, rho2 = thermal_wind.read_density_profile('unequal_depth_step.dens')
    thermal_wind.validate_inputs(
        depth1, depth2, no_motion_depth=280, surface_delta=None)


@nose.tools.raises(ValueError)
def test_no_motion_depth_deeper_than_data():
    depth1, rho1 = thermal_wind.read_density_profile('s109.dens')
    depth2, rho2 = thermal_wind.read_density_profile('s105.dens')
    thermal_wind.validate_inputs(
        depth1, depth1, no_motion_depth=10000, surface_delta=None)


def test_result_in_shallowest_profile_depths():
    depth1, rho1 = thermal_wind.read_density_profile('s109.dens')
    depth2, rho2 = thermal_wind.read_density_profile('s105.dens')
    v_vel_depth = thermal_wind.validate_inputs(
        depth1, depth2, no_motion_depth=280, surface_delta=None)
    np.testing.assert_equal(v_vel_depth, depth1)


def test_depth_result_is_shallowest_profile():
    depth, v_vel = thermal_wind.thermal_wind(
        's109.dens', 's109.dens', no_motion_depth=280)
    expected = np.arange(1, 295, 2)
    np.testing.assert_equal(depth, expected)


def test_v_vel_result_zero_for_degenerate_surface_delta_case():
    depth, v_vel = thermal_wind.thermal_wind(
        's109.dens', 's109.dens', surface_delta=0)
    expected = np.zeros(147)
    np.testing.assert_equal(v_vel, expected)


def test_v_vel_result_zero_for_analytical_no_motion_depth_case():
    del_rho = 2  # kg / m^3
    no_motion_depth = 100  # m
    depth, v_vel = thermal_wind.thermal_wind(
        'constant.dens', 'constant_diff.dens', no_motion_depth=no_motion_depth)
    distance = 36.89e3  # distance between stations in m
    avg_latitude = 26  # average latitude of stations in degrees
    acc_grav = 9.8  # m / s^2
    rho0 = 1024  # kg / m^3
    Omega = 2 * np.pi / (60 * 60 * 24)
    fo = 2 * Omega * np.sin(np.pi * avg_latitude / 180)
    const = acc_grav * del_rho / (rho0 * fo * distance)
    expected = (
        const * (depth - depth ** 2 / (2 * no_motion_depth)
                 - no_motion_depth / 2))
    np.testing.assert_almost_equal(v_vel, expected, 2)

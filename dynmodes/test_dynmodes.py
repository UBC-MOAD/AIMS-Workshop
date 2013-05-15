"""Unit tests for dynmodes module.
"""
import nose.tools
import numpy as np
import dynmodes


def test_build_d2dz2_matrix_n_steps():
    """build_d2dz2_matrix returns number of vertical coordinate grid steps
    """
    depth = np.linspace(0, 1, 3)
    d2dz2, nz, dz = dynmodes.build_d2dz2_matrix(depth)
    nose.tools.assert_equal(nz, 3)


def test_depth2Nsq_Nsq():
    """depth2Nsq returns expected Nsq array
    """
    depth = np.linspace(1, 4, 4)
    density = np.linspace(1028, 1031, 4)
    Nsq_depth, Nsq = dynmodes.density2Nsq(depth, density)
    np.testing.assert_equal(Nsq, np.array([9.8 / 1028] * 4))


def test_depth2Nsq_Nsq_depth():
    """depth2Nsq returns expected Nsq_depth array
    """
    depth = np.linspace(1, 4, 4)
    density = np.linspace(1028, 1031, 4)
    Nsq_depth, Nsq = dynmodes.density2Nsq(depth, density)
    np.testing.assert_equal(Nsq_depth, np.array([0, 1.5, 2.5, 3.5]))


def test_depth2Nsq_Nsq_no_neagative_values():
    """depth2Nsq clamps negative Nsq values to zero
    """
    depth = np.linspace(1, 4, 4)
    density = np.array([1028., 1029., 1030., 1029.])
    Nsq_depth, Nsq = dynmodes.density2Nsq(depth, density)
    np.testing.assert_equal(Nsq, np.array([9.8 / 1028] * 3 + [0]))

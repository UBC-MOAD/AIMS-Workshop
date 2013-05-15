"""Unit tests for dynmodes module.
"""
import nose.tools
import numpy as np
import dynmodes


def test_build_d2dz2_matrix_d2dz2():
    """build_d2dz2_matrix returns 2nd derivative matrix
    """
    depth = np.linspace(0, 1, 4)
    d2dz2, nz, dz = dynmodes.build_d2dz2_matrix(depth)
    expected = np.array((
        (-1,  0,  0,  0),
        (-9, 18, -9,  0),
        ( 0, -9, 18, -9),
        (-1,  0,  0,  0),
    ), dtype=float)
    np.testing.assert_almost_equal(d2dz2, expected)


def test_build_d2dz2_matrix_n_steps():
    """build_d2dz2_matrix returns number of vertical coordinate grid steps
    """
    depth = np.linspace(0, 1, 3)
    d2dz2, nz, dz = dynmodes.build_d2dz2_matrix(depth)
    nose.tools.assert_equal(nz, 3)


def test_build_d2dz2_matrix_dz():
    """build_d2dz2_matrix returns vertical coordinate grid spacing
    """
    depth = np.linspace(0, 1, 3)
    d2dz2, nz, dz = dynmodes.build_d2dz2_matrix(depth)
    np.testing.assert_equal(dz, np.array((0.5, 0.5)))


def test_clean_up_modes_transposes_modes():
    """clean_up_modes returns transpose of modes
    """
    modes = np.array((
        (1, 4),
        (2, 5),
        (3, 6),
    ))
    eigv, modes = dynmodes.clean_up_modes(np.array((7, 8)), modes, 2)
    expected = np.array((
        (1, 2, 3),
        (4, 5, 6),
    ), dtype=float)
    np.testing.assert_equal(modes, expected)


def test_clean_up_modes_excludes_imaginary_eigens():
    """clean_up_modes excludes eigenvalues & modes for imaginary eigenvalues
    """
    eigv = np.array((7j, 8))
    modes = np.array((
        (1, 4),
        (2, 5),
        (3, 6),
    ))
    eigv, modes = dynmodes.clean_up_modes(eigv, modes, 2)
    expected = np.array((
        (4, 5, 6),
    ), dtype=float)
    np.testing.assert_equal(modes, expected)
    np.testing.assert_equal(eigv, np.array((8,)))


def test_clean_up_modes_excludes_tiny_and_negative_eigens():
    """clean_up_modes excludes eigenvalues & modes for near-zero & -ve eigvs
    """
    eigv = np.array((-7, 8, 9e-12))
    modes = np.array((
        (1, 4, 7),
        (2, 5, 8),
        (3, 6, 9),
    ))
    eigv, modes = dynmodes.clean_up_modes(eigv, modes, 3)
    expected = np.array((
        (4, 5, 6),
    ), dtype=float)
    np.testing.assert_equal(modes, expected)
    np.testing.assert_equal(eigv, np.array((8,)))


def test_clean_up_modes_sorts_eigens_by_ascending_eigenvalues():
    """clean_up_modes sorts eigenvalues & modes by ascending eigvalues
    """
    eigv = np.array((9, 7, 8))
    modes = np.array((
        (1, 4, 7),
        (2, 5, 8),
        (3, 6, 9),
    ))
    eigv, modes = dynmodes.clean_up_modes(eigv, modes, 3)
    expected = np.array((
        (4, 5, 6),
        (7, 8, 9),
        (1, 2, 3),
    ), dtype=float)
    np.testing.assert_equal(modes, expected)
    np.testing.assert_equal(eigv, np.array((7, 8, 9)))


def test_clean_up_modes_returns_only_real_parts_of_modes():
    """clean_up_modes returns only real parts of modes
    """
    modes = np.array((
        (1 + 0j, 4),
        (2,      5),
        (3,      6 + 0j),
    ))
    eigv, modes = dynmodes.clean_up_modes(np.array((7, 8)), modes, 2)
    expected = np.array((
        (1, 2, 3),
        (4, 5, 6),
    ), dtype=float)
    np.testing.assert_equal(modes, expected)


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

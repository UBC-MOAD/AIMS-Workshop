"""Ocean vertical dynamic modes exercise code for Mathematical Modelling of Geophysical
Fluids MPE2013 Workshop at African Institute for Mathematical Sciences.
"""
import matplotlib.pyplot as plt
import numpy as np
import scipy.linalg


def dynmodes(Nsq, depth, nmodes):
    """Calculate the 1st nmodes ocean dynamic vertical modes
    given a profile of Brunt-Vaisala (buoyancy) frequencies squared.

    Based on http://woodshole.er.usgs.gov/operations/sea-mat/klinck-html/dynmodes.html
    by John Klinck, 1999.

    :arg Nsq: Brunt-Vaisala (buoyancy) frequencies squared in [1/s^2]
    :type Nsq: :class:`numpy.ndarray`

    :arg depth: Depths in [m]
    :type depth: :class:`numpy.ndarray`

    :arg nmodes: Number of modes to calculate
    :type nmodes: int

    :returns: :obj:`(wmodes, pmodes, rmodes, ce)` vertical velocity modes,
              horizontal velocity modes, vertical density modes, modal speeds
    :rtype: tuple of :class:`numpy.ndarray`
    """
    # 2nd derivative matrix plus boundary conditions
    d2dz2, nz, dz = build_d2dz2_matrix(depth)
    # N-squared diagonal matrix
    Nsq_mat = np.identity(nz) * Nsq
    # Solve generalized eigenvalue problem for eigenvalues and vertical
    # velocity modes
    eigenvalues, wmodes = scipy.linalg.eig(d2dz2, Nsq_mat)
    eigenvalues, wmodes = clean_up_modes(eigenvalues, wmodes, nmodes)
    # Vertical density modes
    rmodes = wmodes * Nsq
    # Modal speeds
    ce = np.real(1 / np.sqrt(eigenvalues))
    # Horizontal velocity modes
    pmodes = np.zeros_like(wmodes)
    # 1st derivative of vertical modes
    pr = np.diff(wmodes) / dz
    # Linear interpolation on to the vertical coordinate grid
    pmodes[:, 1:nz - 1] = (pr[:, 1:nz - 1] + pr[:, :nz - 2]) / 2
    pmodes[:, 0] = pr[:, 0]
    pmodes[:, nz - 1] = pr[:, nz - 2]
    return wmodes, pmodes, rmodes, ce


def build_d2dz2_matrix(depth):
    """Build the matrix that discretizes the 2nd derivative
    over the vertical coordinate, and applies the boundary conditions.

    :arg depth: Depths in [m]
    :type depth: :class:`numpy.ndarray`

    :returns: :obj:`(d2dz2, nz, dz)` 2nd derivative matrix
              (:class:`numpy.ndarray`),
              number of vertical coordinate grid steps,
              and array of vertical coordinate grid point spacings
              (:class:`numpy.ndarray`)
    :rtype: tuple
    """
    zed = -depth
    # Number and size (in [m]) of vertical coordinate grid steps
    nz = zed.size
    dz = zed[:nz - 1] - zed[1:]
    # Grid step midpoints
    z_mid = zed[:nz - 1] - 0.5 * dz[:nz - 1]
    # Size of steps between grid step midpoints
    dz_mid = np.zeros_like(zed)
    dz_mid[1:dz_mid.size - 1] = z_mid[:nz - 2] - z_mid[1:nz - 1]
    dz_mid[0] = dz_mid[1]
    dz_mid[dz_mid.size - 1] = dz_mid[dz_mid.size - 2]
    # 2nd derivative matrix matrix
    d2dz2 = np.zeros((nz, nz))
    # Elements on the diagonal
    diag_index = np.arange(1, nz - 1)
    d2dz2[diag_index, diag_index] = (
        1 / (dz[:nz - 2] * dz_mid[1:nz - 1]) + 1 / (dz[1:nz - 1] * dz_mid[1:nz - 1]))
    # Elements on the super-diagonal
    diag_index_p1 = np.arange(2, nz)
    d2dz2[diag_index, diag_index_p1] = -1 / (dz[:nz - 2] * dz_mid[1:nz - 1])
    # Elements on the sub-diagonal
    diag_index_m1 = np.arange(nz - 2)
    d2dz2[diag_index, diag_index_m1] = -1 / (dz[1:nz - 1] * dz_mid[1:nz - 1])
    # Boundary conditions
    d2dz2[0, 0] = d2dz2[nz - 1, 0] = -1
    return d2dz2, nz, dz


def clean_up_modes(eigenvalues, wmodes, nmodes):
    """Exclude complex-valued and near-zero/negative eigenvalues and their modes.
    Sort the eigenvalues and mode by increasing eigenvalue magnitude,
    truncate the results to the number of modes that were requested,
    and convert the modes from complex to real numbers.

    :arg eigenvalues: Eigenvalues
    :type eigenvalues: :class:`numpy.ndarray`

    :arg wmodes: Modes
    :type wmodes: :class:`numpy.ndarray`

    :arg nmodes: Number of modes requested
    :type nmodes: int

    :returns: :obj:`(eigenvalues, wmodes)`
    :rtype: tuple of :class:`numpy.ndarray`
    """
    # Transpose modes to that they can be handled as an array of vectors
    wmodes = wmodes.transpose()
    # Filter out complex-values and small/negative eigenvalues
    # and corresponding modes
    mask = np.imag(eigenvalues) == 0
    eigenvalues = eigenvalues[mask]
    wmodes = wmodes[mask]
    mask = eigenvalues >= 1e-10
    eigenvalues = eigenvalues[mask]
    wmodes = wmodes[mask]
    # Sort eigenvalues and modes and truncate to number of modes requests
    index = np.argsort(eigenvalues)
    eigenvalues = eigenvalues[index[:nmodes]]
    wmodes = wmodes[index[:nmodes]]
    return eigenvalues, wmodes.real


def plot_modes(Nsq, depth, nmodes, wmodes, pmodes, rmodes):
    """Plot Brunt-Vaisala (buoyancy) frequency profile and 3 sets of modes
    (vertical velocity, horizontal velocity, and vertical density) in 4 panes.

    :arg Nsq: Brunt-Vaisala (buoyancy) frequencies squared in [1/s^2]
    :type Nsq: :class:`numpy.ndarray`

    :arg depth: Depths in [m]
    :type depth: :class:`numpy.ndarray`

    :arg wmodes: Vertical velocity modes
    :type wmodes: :class:`numpy.ndarray`

    :arg pmodes: Horizontal velocity modes
    :type pmodes: :class:`numpy.ndarray`

    :arg rmodes: Vertical density modes
    :type rmodes: :class:`numpy.ndarray`

    :arg nmodes: Number of modes to calculate
    :type nmodes: int
    """
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(2, 2, 1)
    # Nsq
    ax.plot(Nsq, -depth)
    ax.ticklabel_format(style='sci', scilimits=(2, 2), axis='x')
    ax.set_ylabel('z')
    ax.set_xlabel('N^2')
    # modes
    mode_sets = [
        # (values, subplot number, x-axis title)
        (wmodes, 2, 'wmodes'),
        (pmodes, 3, 'pmodes'),
        (rmodes, 4, 'rmodes'),
    ]
    for mode_set in mode_sets:
        modes, subplot, title = mode_set
        ax = fig.add_subplot(2, 2, subplot)
        for i in xrange(nmodes):
            ax.plot(modes[i], -depth, label='mode {}'.format(i + 1))
        ax.ticklabel_format(style='sci', scilimits=(3, 3), axis='x')
        ax.set_ylabel('z')
        ax.set_xlabel(title)
        ax.legend(loc='best')


def read_density_profile(filename):
    """Return depth and density arrays read from filename.

    :arg filename: Name of density profile file.
    :type filename: string

    :returns: :obj:`(depth, density)` depths, densities
    :rtype: tuple of :class:`numpy.ndarray`
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


def density2Nsq(depth, density, rho0=1028):
    """Return the Brunt-Vaisala (buoyancy) frequency (Nsq) profile
    corresponding to the given density profile.
    The surface Nsq value is set to the value of the 1st calculated value
    below the surface.
    Also return the depths for which the Brunt-Vaisala (buoyancy) frequencies squared
    were calculated.

    :arg depth: Depths in [m]
    :type depth: :class:`numpy.ndarray`

    :arg density: Densities in [kg/m^3]
    :type density: :class:`numpy.ndarray`

    :arg rho0: Reference density in [kg/m^3]; defaults to 1028
    :type rho0: number

    :returns: :obj:`(Nsq_depth, Nsq)` depths for which the Brunt-Vaisala
              (buoyancy) frequencies squared were calculated,
              Brunt-Vaisala (buoyancy) frequencies squared
    :rtype: tuple of :class:`numpy.ndarray`
    """
    grav_acc = 9.8  # m / s^2
    Nsq = np.zeros_like(density)
    Nsq[1:] = np.diff(density) * grav_acc / (np.diff(depth) * rho0)
    Nsq[0] = Nsq[1]
    Nsq[Nsq < 0] = 0
    Nsq_depth = np.zeros_like(depth)
    Nsq_depth[1:] = (depth[:depth.size - 1] + depth[1:]) / 2
    return Nsq_depth, Nsq

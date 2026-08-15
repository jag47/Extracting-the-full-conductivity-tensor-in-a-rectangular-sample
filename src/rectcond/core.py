"""Extraction of the full 2D conductivity tensor in a rectangular sample.

Implements the procedure of J. Gelfond and O. Vafek,
"Extracting the full conductivity tensor in a rectangular sample",
arXiv:2608.06643 (2026). Equation numbers in comments refer to that paper.

This module is a direct transcription of the original release script
("extract conductivity tensor - python.py" in this repository), wrapped
in a function so it can be called on measured data.
"""

from dataclasses import dataclass

import numpy as np
import scipy.optimize as opt
import scipy.special as sp
from mpmath import mp

__all__ = ["ConductivityTensor", "extract_tensor"]


@dataclass
class ConductivityTensor:
    """Result of the conductivity-tensor extraction.

    Conductivities are sheet conductances in units of 1/[R], where [R] is
    the unit of the input resistances. ``alpha`` is in radians.
    """

    sigma_plus: float
    """Larger principal conductivity (sigma_+ >= sigma_-)."""
    sigma_minus: float
    """Smaller principal conductivity."""
    sigma_hall: float
    """Hall conductivity sigma_H."""
    alpha: float
    """Principal-axes angle in [0, pi), counterclockwise from the bottom
    (d1) edge of the sample."""
    geometric_mean: float
    """sqrt(sigma_+ * sigma_-)."""
    anisotropy_ratio: float
    """sqrt(sigma_- / sigma_+)."""
    r: float
    """Hypergeometric parameter r of the Schwarz-Christoffel map."""
    a: float
    """Hypergeometric parameter a; a = 1/2 corresponds to principal axes
    aligned with the sample edges."""
    z5: float
    """Image of the bottom-edge midpoint in the complex z-plane."""


def _k(a, r):
    # K_a(r) up to the pi/2 prefactor, Eq. (13)
    return sp.hyp2f1(a, 1 - a, 1, r**2)


def extract_tensor(d1, d2, R1, R2, R3, R5):
    """Extract the full conductivity tensor from four resistance measurements.

    The sample is a rectangle of side d1 along x (horizontal) and d2 along
    y (vertical), bottom-left corner at the origin. Each resistance is a
    four-terminal measurement (Phi_A - Phi_B)/I with current I from source
    S to drain D, in the configurations of Fig. 1f-h of the paper:

    ============= ============ =========== ================ ==========
    configuration source S     drain D     probe A (V+)     probe B (V-)
    ============= ============ =========== ================ ==========
    1             bottom-left  bottom-right top-right       top-left
    2             bottom-right top-right   bottom-left      top-left
    3             bottom-left  top-right   bottom-right     top-left
    5             bottom-right top-right   bottom-edge mid  top-left
    ============= ============ =========== ================ ==========

    R1 is negated so that it is positive, R1 = -(Phi_A - Phi_B)/I; the
    others keep their sign (R3 in particular may be negative, and its sign
    carries the sign of the Hall response). Consistency requires R2 < R5.
    A fourth vertex configuration (S bottom-right, D top-left, A
    bottom-left, B top-right) is not needed but provides the cross-check
    R4 = -2*R1 + 2*R2 - R3, Eq. (39).

    Parameters
    ----------
    d1, d2 : float
        Side lengths of the sample (any common unit; only the ratio enters).
    R1, R2, R3 : float
        Vertex-configuration resistances (configurations 1-3).
    R5 : float
        Midpoint-configuration resistance (configuration 5).

    Returns
    -------
    ConductivityTensor
        Principal conductivities sigma_+ >= sigma_-, Hall conductivity,
        principal-axes angle alpha, and the intermediate parameters
        (r, a, z5) for diagnostics.
    """
    # find r, hypergeometric parameter: log form of 1 - r^2 = r^(2 R1/R2), Eq. (31)
    def find_r(x):
        return np.log(1 / x**2) / np.log(1 - x**2) + R2 / R1

    sol = opt.root(find_r, 0.5)
    r = sol.x[0]

    # Hall conductivity, Eq. (40)
    sigma_hall = (np.log(1 - r**2) * (R3 * np.log(1 - r**2) + R1 * np.log(1 / r**2 - 1)) /
                  ((R1 * np.pi)**2 + (R1 * np.log(1 / r**2 - 1) + R3 * np.log(1 - r**2))**2))

    # geometric mean sqrt(sigma_+ sigma_-), Eq. (41)
    geometric_mean = (-R1 * np.log(1 - r**2) * np.pi /
                      ((R1 * np.pi)**2 + (R1 * np.log(1 / r**2 - 1) + R3 * np.log(1 - r**2))**2))

    # image z5 of the bottom-edge midpoint in the z-plane, Eq. (42)
    z5 = ((1 / r**2 - (1 - r**2)**(-R5 / R1)) /
          (1 - (1 - r**2)**(-R5 / R1)))

    # find a, hypergeometric parameter: midpoint constraint, Eq. (44)
    def find_a(x):
        return (np.pi / 2 * mp.hyp2f1(x, 1 - x, 1, r**2) * (1 - x)
                - mp.appellf1(1 - x, 1 - x, x, 2 - x, z5, z5 * r**2) * z5**(1 - x) * mp.sin(np.pi * x))

    a = np.float64(mp.findroot(find_a, mp.mpf(0.25)))

    # principal-axes angle alpha, Eq. (45); branch chosen so that
    # a < 1/2 -> alpha in (0, pi/2) and a > 1/2 -> alpha in (pi/2, pi)
    alpha = 1 / 2 * np.arctan(2 * d1 * d2 * _k(a, r) * _k(a, np.sqrt(1 - r**2)) /
                              (d1**2 * _k(a, np.sqrt(1 - r**2))**2 - d2**2 * _k(a, r)**2) *
                              np.cos(np.pi * a))
    if a > 1 / 2:
        while alpha < np.pi / 2:
            alpha += np.pi / 2
    else:
        while alpha < 0:
            alpha += np.pi / 2

    # anisotropy ratio sqrt(sigma_- / sigma_+), Eq. (47)
    anisotropy_ratio = ((-1 + np.sqrt(1 + np.tan(a * np.pi)**2 * np.sin(2 * alpha)**2)) /
                        (np.tan(a * np.pi) * np.sin(2 * alpha)))

    sigma_minus = geometric_mean * anisotropy_ratio
    sigma_plus = geometric_mean / anisotropy_ratio

    return ConductivityTensor(
        sigma_plus=float(sigma_plus),
        sigma_minus=float(sigma_minus),
        sigma_hall=float(sigma_hall),
        alpha=float(alpha),
        geometric_mean=float(geometric_mean),
        anisotropy_ratio=float(anisotropy_ratio),
        r=float(r),
        a=float(a),
        z5=float(z5),
    )

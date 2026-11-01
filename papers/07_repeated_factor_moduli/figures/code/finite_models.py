"""Explicit finite forms and the polynomial curve; no HE solver."""
from functools import lru_cache
import numpy as np
import sympy as sp


MASS_VECTOR = np.array([1., 1., 0., 0., 1., 1.])
MASS = np.eye(6) + np.outer(MASS_VECTOR, MASS_VECTOR)
WHITEN = np.linalg.inv(np.linalg.cholesky(MASS).T)


def blocks(coordinates):
    z1, m11, m12, m21, m22, z3 = coordinates
    return z1, np.array([[m11, m12], [m21, m22]]), z3, -z1-m11-m22-z3


def balanced_rows(epsilon):
    e = epsilon**2
    z = (4*e-2+np.sqrt(64+44*e+16*e**2))/10
    rows = []
    for coordinate in np.eye(6):
        z1, matrix, z3, z4 = blocks(coordinate)
        row = np.r_[
            np.sqrt(1+e)*(np.array([z1, 0.])-matrix[0]),
            np.sqrt(e)*(np.array([0., z1])-matrix[1]),
            np.sqrt((1+e)/z)*(z1-z3),
            matrix[:, 0]-np.array([z4, 0.]),
            np.sqrt(z)*(z3-z4),
        ]
        rows.append(row)
    return np.array(rows).T


def eigenvalues(epsilon):
    # SVD of the weighted differential avoids subtractive cancellation
    # in the smallest eigenvalue of its normal matrix.
    singular = np.linalg.svd(balanced_rows(epsilon) @ WHITEN,
                             compute_uv=False)
    return np.sort(singular**2)


@lru_cache(maxsize=2)
def curve_functions(rotated=True):
    epsilon, t = sp.symbols("epsilon t", positive=True, real=True)
    e = epsilon**2
    z = (4*e-2+sp.sqrt(64+44*e+16*e**2))/10
    temp = t**2/(18*sp.pi)
    aa, bb = e+1+z, -2*e+7*(1+z)
    beta = 1/(1-temp*bb)
    p = sp.Matrix([5*temp*beta*aa, beta-1])
    alpha = 1+p[0]
    ll = sp.sqrt(beta/alpha)/sp.sqrt(2)*sp.Matrix(
        [[sp.sqrt(1+e), epsilon], [sp.sqrt(1+e), -epsilon]])
    pol = sp.Matrix([sp.cos(epsilon), 0, sp.sin(epsilon)])
    if not rotated:
        pol = sp.Matrix([1, 0, 0])
    qq = sp.Matrix([0, 1, 0])
    row = t*(pol*ll[0, :] + qq*ll[1, :])
    scalar = t*sp.sqrt((1+e)/(2*z))*(pol+qq)
    col = sp.Matrix([t, 0])
    final = t*sp.sqrt(beta*z/alpha)
    coeff = sp.Matrix(list(row)+list(scalar)+list(col)+[final])
    return sp.lambdify((t, epsilon),
        (coeff, coeff.diff(epsilon), p, p.diff(epsilon), alpha, beta),
        modules="numpy", cse=True)


def curve_data(t, epsilon, rotated=True):
    coeff, vel, p, dp, alpha, beta = curve_functions(rotated)(t, epsilon)
    coeff, vel = np.asarray(coeff).ravel(), np.asarray(vel).ravel()
    row, scalar, col, final = (coeff[:6].reshape(3, 2), coeff[6:9],
                               coeff[9:11], coeff[11])
    orbit = []
    for coordinate in np.eye(6):
        z1, matrix, z3, z4 = blocks(coordinate)
        orbit.append(np.r_[(z1*row-row@matrix).ravel(),
                            (z1-z3)*scalar, matrix@col-col*z4,
                            (z3-z4)*final])
    orbit = np.array(orbit).T
    weights = np.sqrt(np.r_[np.full(6, alpha), np.full(5, beta), alpha])
    weighted_orbit = weights[:, None]*orbit
    weighted_vel = weights*vel
    gauge = np.linalg.lstsq(weighted_orbit, weighted_vel, rcond=1e-13)[0]
    remainder = weighted_vel-weighted_orbit@gauge
    finite = float(remainder@remainder)
    parameter = float(np.asarray(dp).ravel()@np.asarray(dp).ravel())
    lower = t*t*beta*(1+2*epsilon**2)/2 if rotated else 0.
    return finite, parameter, lower, np.asarray(p).ravel(), coeff

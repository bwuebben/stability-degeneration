"""Exact finite controls and numerical balancing for the rank-five example.

Run with: uv run --no-project --with sympy==1.14.0 --with numpy==2.5.4
  --with scipy==1.18.1 python figures/code/verify_star_model.py
Finite checks only: this script does not solve or certify the HE PDE.
"""
from pathlib import Path
import json
import sys

import numpy as np
import scipy
from scipy.linalg import null_space
from scipy.optimize import root
import sympy as sp


B_EXACT = [
    sp.Matrix([[sp.Rational(3, 2), 0], [0, sp.Rational(1, 2)]]),
    sp.Matrix([[1, sp.Rational(1, 2)], [sp.Rational(1, 2), 1]]),
    sp.Matrix([[sp.Rational(1, 2), -sp.Rational(1, 2)],
               [-sp.Rational(1, 2), sp.Rational(3, 2)]]),
]
Z_EXACT = [
    sp.Matrix([[sp.sqrt(6)/2, 0, 0], [0, sp.sqrt(2)/2, 0]]),
    sp.Matrix([[1, 0, 0], [sp.Rational(1, 2), sp.sqrt(3)/2, 0]]),
    sp.Matrix([[sp.sqrt(2)/2, 0, 0], [-sp.sqrt(2)/2, 1, 0]]),
]


def matrix_operator(covs):
    """Full complex commutator form on M_2(C) plus three scalar vertices."""
    total = sum(covs, sp.zeros(2))
    cols = []
    for k in range(7):
        x = sp.zeros(7, 1)
        x[k] = 1
        target = sp.Matrix(2, 2, list(x[:4]))
        sources = x[4:]
        target_out = target * total - sum(
            (sources[j] * covs[j] for j in range(3)), sp.zeros(2))
        source_out = [sources[j] * sp.trace(covs[j])
                      - sp.trace(covs[j] * target) for j in range(3)]
        cols.append(sp.Matrix(list(target_out) + source_out))
    return sp.Matrix.hstack(*cols)


def exact_checks():
    assert sum(B_EXACT, sp.zeros(2)) == 3 * sp.eye(2)
    for b, z in zip(B_EXACT, Z_EXACT):
        assert sp.simplify(z * z.T - b) == sp.zeros(2)
        assert sp.trace(b) == 2 and b.det() > 0 and b[0, 0] > 0
        assert z.rank() == 2
    comm = B_EXACT[0] * B_EXACT[1] - B_EXACT[1] * B_EXACT[0]
    assert comm != sp.zeros(2)
    lap = matrix_operator(B_EXACT)
    ident = sp.Matrix([1, 0, 0, 1, 1, 1, 1])
    assert lap == lap.T and lap * ident == sp.zeros(7, 1)
    assert lap.rank() == 6
    embed = sp.eye(7)[:, :6]
    for k in range(6):
        embed[6, k] = -ident[k]
    mass = embed.T * embed
    gram = embed.T * lap * embed
    lam = sp.Symbol('lambda')
    poly = sp.factor((gram - lam * mass).det() / mass.det())
    # Every leading principal minor is positive: exact positivity on trace-free space.
    minors = [gram[:k, :k].det() for k in range(1, 7)]
    assert all(x > 0 for x in minors)
    # Scalar-vertex compression has only three modes and misses half the full space.
    scalar_embed = sp.Matrix([[1,0,0,0], [0,0,0,0], [0,0,0,0],
                              [1,0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,1]])
    scalar_mass = scalar_embed.T * scalar_embed
    scalar_gram = scalar_embed.T * lap * scalar_embed
    scalar_poly = sp.factor((scalar_gram-lam*scalar_mass).det()/scalar_mass.det())
    # Central balance does preserve adjoint (the matrices here are real).
    transpose = sp.eye(7)
    transpose[1,1] = transpose[2,2] = 0
    transpose[1,2] = transpose[2,1] = 1
    assert lap*transpose == transpose*lap
    # Off balance, a Hermitian-only kernel test can miss a complex nilpotent.
    jordan = sp.Matrix([[0, 1], [0, 0]])
    arrows = [sp.eye(2), jordan]
    raw_cols = []
    for k in range(8):
        unit = sp.eye(8)[:,k]
        x = sp.Matrix(2,2,list(unit[:4]))
        y = sp.Matrix(2,2,list(unit[4:]))
        raw_cols.append(sp.Matrix([entry for z in arrows for entry in x*z-z*y]))
    raw = sp.Matrix.hstack(*raw_cols)
    jordan_trace = sp.Matrix([1,0,0,1,1,0,0,1])
    jordan_embed = sp.eye(8)[:,:7]
    for k in range(7):
        jordan_embed[7,k] = -jordan_trace[k]
    complex_rank = (raw*jordan_embed).rank()
    pauli = [sp.Matrix([[0,1],[1,0]]),sp.Matrix([[0,-sp.I],[sp.I,0]]),
             sp.diag(1,-1)]
    hermitian_cols = [sp.Matrix(list(sp.eye(2))+list(-sp.eye(2)))]
    hermitian_cols += [sp.Matrix(list(p)+list(sp.zeros(2))) for p in pauli]
    hermitian_cols += [sp.Matrix(list(sp.zeros(2))+list(p)) for p in pauli]
    hermitian_derivative = raw*sp.Matrix.hstack(*hermitian_cols)
    hermitian_gram = (hermitian_derivative.H*hermitian_derivative).applyfunc(sp.re)
    assert complex_rank == 6 and hermitian_gram.rank() == 7
    # A connected identity-arrow representation has a four-dimensional centralizer.
    identity_arrow = raw[:4,:]
    assert 8-identity_arrow.rank() == 4
    moment_target = sum((z*z.H for z in arrows),sp.zeros(2))
    moment_source = -sum((z.H*z for z in arrows),sp.zeros(2))
    assert moment_target == sp.diag(2,1)
    assert moment_source == sp.diag(-1,-2)
    return {'covariances': [str(b) for b in B_EXACT],
            'determinants': [str(b.det()) for b in B_EXACT],
            'noncommuting_covariance_commutator': str(comm),
            'full_trace_free_dimension': 6,
            'centralizer_dimension': 1,
            'identity_arrow_centralizer_dimension': 4,
            'trace_free_characteristic_polynomial': str(poly),
            'positive_gram_leading_minors': [str(x) for x in minors],
            'scalar_compression_characteristic_polynomial_including_zero': str(scalar_poly),
            'central_balance_adjoint_invariance': True,
            'off_balance_jordan_control': {
                'complex_trace_free_dimension': 7,
                'complex_trace_free_kernel_dimension': 7-complex_rank,
                'hermitian_real_kernel_dimension': 7-hermitian_gram.rank(),
                'hermitian_gram_determinant': str(hermitian_gram.det()),
                'target_moment': str(moment_target),'source_moment': str(moment_source)}}


B_NUM = [np.array(b, dtype=float) for b in B_EXACT]
TRACE_FREE = null_space(np.array([[1., 0., 0., 1., 1., 1., 1.]]))


def positive_gauge(y):
    x, z = y
    radius = np.hypot(x, z)
    s = 0.5 if radius < 1e-14 else np.sinh(radius/2)/radius
    return np.cosh(radius/2)*np.eye(2) + s*np.array([[x, z], [z, -x]])


def balanced_covariances(weights):
    weights = np.asarray(weights, dtype=float)
    normalized = weights/weights.sum()
    def covariances(y):
        gauge = positive_gauge(y)
        mats = [gauge @ b @ gauge for b in B_NUM]
        return [a*b/np.trace(b) for a,b in zip(normalized,mats)]
    def equations(y):
        total = sum(covariances(y))
        return [total[0,0]-total[1,1], 2*total[0,1]]
    solve = root(equations, np.zeros(2), tol=1e-11)
    error = float(np.linalg.norm(equations(solve.x)))
    assert error < 1e-10, (solve.message,error)
    covs = [weights.sum()*b for b in covariances(solve.x)]
    return covs, solve.x, error


def full_spectrum(covs):
    """Use the unsquared complex commutator map to retain tiny modes."""
    arrows = [np.linalg.cholesky(c) for c in covs]
    columns=[]
    for k in range(7):
        unit=np.eye(7)[:,k]
        target=unit[:4].reshape(2,2)
        columns.append(np.concatenate([
            (target@z-unit[4+j]*z).reshape(-1) for j,z in enumerate(arrows)]))
    deriv=np.column_stack(columns)@TRACE_FREE
    return np.sort(np.linalg.svd(deriv,compute_uv=False)**2)


def numerical_checks():
    records=[]
    for delta in [1., .1, .01, 1e-4, 1e-6, 1e-8, 1e-10]:
        covs,gauge,error=balanced_covariances([delta,1.,1.])
        vals=full_spectrum(covs)
        assert np.all(vals>0)
        assert np.linalg.norm(sum(covs)-(delta+2)/2*np.eye(2)) < 1e-9
        if delta<=1e-4:
            assert abs(vals[0]/delta-1.25)<max(2*delta,2e-6)
        records.append({'weak_demand':delta,'gauge_coordinates':gauge.tolist(),
                        'normalized_balance_residual':error,
                        'eigenvalues':vals.tolist(),
                        'weak_coefficient':float(vals[0]/delta),
                        'next_eigenvalue':float(vals[1]),
                        'covariance_commutator_norm':float(np.linalg.norm(covs[0]@covs[1]-covs[1]@covs[0]))})
    return records


if __name__=='__main__':
    result={'scope':'Finite exact model checks and floating-point balancing; no PDE verification.',
            'versions':{'python':sys.version.split()[0],'sympy':sp.__version__,
                        'numpy':np.__version__,'scipy':scipy.__version__},
            'exact':exact_checks(),'numerical':numerical_checks()}
    output=Path(__file__).resolve().parents[1]/'data'/'star_model_checks.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'output':str(output),'exact':result['exact'],
                      'weak_coefficients':[r['weak_coefficient'] for r in result['numerical']]},indent=2))

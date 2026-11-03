"""Supplementary numerical checks for the general-cancellation revision.

This independently differentiates the exact scalar mass expression at high
precision, checks the analytic Hessian used in the figure, and checks the
printed profile derivative constants. It is not a proof of uniform bounds,
flow asymptotics, or Holder failure.
"""
from pathlib import Path
import hashlib
import importlib.util
import json
import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
mp.mp.dps = 100


def exact_mass(epsilon, lam):
    p = mp.exp(-lam**2)
    rho = -mp.expm1(-lam**2)/lam**2
    y = rho**2*epsilon**2/(2+epsilon**2)
    # Sum of the singular values of each 2 by 2 source block is
    # sqrt(tr(Z*Z) + 2|det Z|); the other two blocks are scalar.
    source = mp.sqrt(1+5*p*p+2*y+2*mp.sqrt(5*p*p+2*p*y+y*y))
    return 2*source+2*mp.sqrt((3*rho-2*p)**2+p*p)+2*mp.sqrt(1+9*rho*rho)


def main():
    figure_path = ROOT/'figures/code/fig07_mass_profile.py'
    spec = importlib.util.spec_from_file_location('figure_profile', figure_path)
    figure = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(figure)
    xi, _, _, rows = figure.compute_data()
    errors = []
    for row in rows:
        lam = mp.mpf(str(row['lam']))
        p = mp.exp(-lam**2)
        rho = -mp.expm1(-lam**2)/lam**2
        scale = mp.sqrt(2*p)/rho
        for index in (0, 125, 250, 375, 500):
            epsilon = scale*mp.mpf(str(xi[index]))
            hessian = mp.diff(lambda e: exact_mass(e, lam), epsilon, 2)/rho**2
            errors.append(abs(float(hessian)-row['hessian'][index]))
    assert max(errors) < 1e-12, max(errors)

    def profile(x):
        return x*x+mp.sqrt((1+x*x)**2+4)-mp.sqrt(5)
    expected = [(0, 2, 2+2/mp.sqrt(5)), (1, 2, 2+3/mp.sqrt(2)),
                (1, 3, 3*mp.sqrt(2)/4)]
    for point, order, value in expected:
        assert abs(mp.diff(profile, point, order)-value) < mp.mpf('1e-90')
    delta = 3/mp.sqrt(2)-2/mp.sqrt(5)
    assert delta > 0

    # Independent differentiation of the unscaled exact mass, including the
    # very small scale. These samples check transcription, not uniform limits.
    ratios = []
    for lam in map(mp.mpf, (4, 6, 8)):
        p = mp.exp(-lam**2)
        rho = -mp.expm1(-lam**2)/lam**2
        scale = mp.sqrt(2*p)/rho
        d0 = mp.diff(lambda e: exact_mass(e, lam), 0, 2)
        d1 = mp.diff(lambda e: exact_mass(e, lam), scale, 2)
        ratio = (d1-d0)/(rho*rho*delta)
        ratios.append({'lambda':str(lam), 'two_point_ratio':mp.nstr(ratio, 30)})
    assert abs(mp.mpf(ratios[-1]['two_point_ratio'])-1) < mp.mpf('1e-20')
    report = {
        'scope':'Independent finite/high-precision transcription checks; not proofs of theorems.',
        'precision_decimal_digits':100,
        'figure_hessian_comparisons':len(errors),
        'maximum_absolute_error':max(errors),
        'profile_constants':'Psecond(0)=2+2/sqrt(5); Psecond(1)=2+3/sqrt(2); Pthird(1)=3sqrt(2)/4',
        'positive_hessian_difference':mp.nstr(delta, 40),
        'two_point_checks':ratios,
        'figure_script_sha256':hashlib.sha256(figure_path.read_bytes()).hexdigest(),
        'validation_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    (ROOT/'validation/revision-report.json').write_text(json.dumps(report, indent=2)+'\n')
    print('PASS: 15 exact-mass/figure Hessian comparisons, profile constants, and shrinking-scale two-point checks.')


if __name__ == '__main__':
    main()

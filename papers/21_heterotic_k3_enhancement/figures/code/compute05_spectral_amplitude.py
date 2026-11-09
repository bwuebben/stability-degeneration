"""Exact normalization and cancellation checks for the spectral tree formula.

Run from the paper directory: python figures/code/compute05_spectral_amplitude.py
These finite identities do not test analytic uniformity or a KK-tail estimate.
"""

import json
from pathlib import Path

import sympy as sp


def zero(value):
    assert sp.simplify(value) == 0, value


g10, L, lam, mass2, volume, a, b = sp.symbols(
    'g10 L lam mass2 V0 a b', positive=True)
phi, chi, normal = sp.symbols('phi chi normal', real=True)
buu, buv, bvv, current = sp.symbols('buu buv bvv j', real=True)
s, t, u, scale, shift = sp.symbols('s t u scale shift', real=True)
gstar = g10 / L**2
mass = sp.sqrt(lam) / L
moment = (sp.sqrt(lam) * g10/L * normal
          + (g10/L)**2 / 2 * (buu*phi**2 + 2*buv*phi*chi + bvv*chi**2))
potential = moment**2 / (2*g10**2)
origin = {phi: 0, chi: 0, normal: 0}
zero(sp.diff(potential, normal, 2).subs(origin) - lam/L**2)
zero(sp.diff(potential, normal, phi, chi).subs(origin) - gstar*mass*buv)
zero(sp.diff(potential, normal, phi, phi).subs(origin) - gstar*mass*buu)
zero(sp.diff(potential, phi, phi, chi, chi).subs(origin)
     - gstar**2*(buu*bvv + 2*buv**2))

contact = -gstar**2*(buu*bvv + 2*buv**2)
exchange = gstar**2*mass2*(buv**2/(mass2-s)
    + buu*bvv/(mass2-t) + buv**2/(mass2-u))
subtracted = gstar**2*(buv**2*(s/(mass2-s)+u/(mass2-u))
    + buu*bvv*t/(mass2-t))
zero(contact+exchange-subtracted)
vector = gstar**2*current**2*((t-u)/(mass2-s)+(t-s)/(mass2-u))
low = sp.series((vector+subtracted).subs(
    {s:scale*s, t:scale*t, u:scale*u}), scale, 0, 2).removeO().coeff(scale)
expected = gstar**2/mass2*(3*current**2+buu*bvv-buv**2)*t
zero((low-expected).subs(u, -s-t))

jlight2 = sp.Rational(2, 3)/volume
mu_light = 8*b**2/(3*volume)
Kreal = 3*jlight2/mu_light
Kphysical = g10**2/L**2*Kreal
e2 = gstar**2*jlight2
zero(Kreal-sp.Rational(3, 4)/b**2)
zero(Kphysical-3*e2/(mu_light/L**2))
zero(1/(mass2+shift-s)-1/(mass2-s)
     + shift/((mass2+shift-s)*(mass2-s)))

record = {
    'status':'PASS',
    'sympy_version':sp.__version__,
    'scope':'Exact finite identities; no numerical proof of analytic or KK bounds.',
    'scalar_mass_squared':str(lam/L**2),
    'normal_scalar_cubic':str(gstar*mass*buv),
    'quartic_potential_vertex':str(gstar**2*(buu*bvv+2*buv**2)),
    'scalar_contact_plus_exchange':str(subtracted),
    'low_energy_curvature_coefficient':str(expected/t),
    'real_reference_curvature':str(sp.simplify(Kreal)),
    'physical_reference_curvature':str(sp.simplify(Kphysical)),
    'checks':['Canonical scalar mass and cubic normalization',
              'Quartic potential symmetry factors',
              'Exact contact/exchange cancellation',
              'All three scalar exchange channels',
              'Gauss plus ONeill low-energy coefficient',
              'Real/Dolbeault and physical metric factors',
              'Pole-location resolvent difference']
}
output = Path(__file__).resolve().parents[1] / 'data' / 'spectral_amplitude.json'
output.parent.mkdir(parents=True, exist_ok=True)
output.write_text(json.dumps(record, indent=2)+'\n')
print(json.dumps(record, indent=2))

"""Exact checks of the printed toric intersections and energy coefficients.
Run from any directory. No data files or numerical HE solutions are used.
"""
from fractions import Fraction as F
import json

def pair(a,b):
    return -a[0]*b[0]+a[0]*b[1]+a[1]*b[0]-a[2]*b[2]
L=(2,1,-1);M=(0,2,0);alpha=(4,7,-2);f=(0,1,0)
c1=tuple(a+b for a,b in zip(L,M))
rays=[(1,0),(1,1),(0,1),(-1,0),(1,-1)]
dets=[a[0]*b[1]-a[1]*b[0] for a,b in zip(rays,rays[1:]+rays[:1])]
values=dict(L2=pair(L,L),M2=pair(M,M),LM=pair(L,M),c1_squared=pair(c1,c1),
    c2=len(rays),degree_L=pair(L,alpha),degree_M=pair(M,alpha),
    polarization_squared=pair(alpha,alpha),volume=pair(alpha,alpha)//2,
    degree_L_derivative=pair(L,f),degree_M_derivative=pair(M,f),
    volume_derivative=pair(alpha,f),ray_determinants=dets,
    tracefree_energy_pi2=8*(F(len(rays))-F(pair(c1,c1),4)),
    limit_energy_pi2=8*(F(pair(L,M))-F(pair(c1,c1),4)))
expected={'L2':-1,'M2':0,'LM':4,'c1_squared':7,'c2':5,'degree_L':8,'degree_M':8,
    'polarization_squared':36,'volume':18,'degree_L_derivative':2,'degree_M_derivative':0,
    'volume_derivative':4,'ray_determinants':[1]*5,'tracefree_energy_pi2':26,'limit_energy_pi2':18}
assert values==expected,(values,expected)
# The minimizer of w_L^2+w_M^2 under w_M-w_L>=1 has mean zero.
wL,wM=F(-1,2),F(1,2)
assert wM-wL==1 and wL+wM==0
# Coefficient norm |beta_0|^2=r^-6 and area(S^3_r)=2*pi^2*r^3:
# integrating 2*r^-3 gives rho^-2-R^-2, coefficient pi^2.
assert -2*F(-1,2)==1
print(json.dumps(values,default=str,indent=2))
print('PASS: printed intersections, volume, slopes, energies and two weights.')

"""Finite-model boundary speed and length bounds; no numerical HE solution."""
from pathlib import Path
import csv
import numpy as np
import sympy as s
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

out=Path(__file__).resolve().parents[1]/'output';out.mkdir(exist_ok=True)
r=s.symbols('r',real=True,nonnegative=True)
c=s.sqrt(s.Rational(3,5))
chi=(s.Rational(6,5))/(s.sqrt((s.Rational(8,5)*r**2)**2+s.Rational(12,5)*(1+r**4))+s.Rational(8,5)*r**2)
f=s.sqrt((1-r**2*chi)/(1+2*r**4))
curve=s.Matrix([f*r**2,f,f*r**2,r*s.sqrt(chi),s.sqrt(s.Rational(3,5)-r**2*chi),r*s.sqrt(chi)])
velocity=curve.diff(r)
speed_fn=s.lambdify(r,s.sqrt(sum(v*v for v in velocity)),'numpy',cse=True)
chi_fn=s.lambdify(r,chi,'numpy',cse=True)
rr=np.linspace(0,.6,2401)
speed=np.asarray(speed_fn(rr));cc=float(c)
upper=np.concatenate(([0],np.cumsum((speed[1:]+speed[:-1])*np.diff(rr)/2)))
lower=rr*np.sqrt(2*chi_fn(rr));leading=rr*np.sqrt(2*cc)
assert np.all(upper+1e-12>=lower)
assert abs(speed[0]**2-2*cc)<1e-12
with (out/'boundary_geometry.csv').open('w',newline='') as file:
 writer=csv.writer(file);writer.writerow(['r','speed_squared_over_l','finite_length_upper_over_sqrt_l','finite_distance_lower_over_sqrt_l','leading_distance_over_sqrt_l'])
 writer.writerows(zip(rr,speed**2,upper,lower,leading))
plt.rcParams.update({'font.family':'serif','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
fig,axs=plt.subplots(1,2,figsize=(7.4,2.9),constrained_layout=True)
axs[0].plot(rr,speed**2/(2*cc),color='#245c9f',lw=1.8)
axs[0].axhline(1,color='#777777',ls=':',lw=1)
axs[0].set(xlabel=r'$r=M^{-1/2}$',ylabel=r"$\|z'(r)\|^2/(2l\sqrt{3/5})$",title='Finite horizontal speed',xlim=(0,.6))
axs[1].plot(rr[1:],upper[1:]/leading[1:],label='Path length upper bound',color='#245c9f',lw=1.8)
axs[1].plot(rr[1:],lower[1:]/leading[1:],label='Charged norm lower bound',color='#b2672f',lw=1.8)
axs[1].axhline(1,color='#777777',ls=':',lw=1)
axs[1].set(xlabel=r'$r=M^{-1/2}$',ylabel=r'Bound divided by $\sqrt{2l\sqrt{3/5}}\,r$',title='Matching distance bounds',xlim=(0,.6))
axs[1].legend(loc='lower left',frameon=False,fontsize=8)
for ax in axs:ax.grid(alpha=.18)
fig.savefig(out/'boundary_geometry.pdf');fig.savefig(out/'boundary_geometry.png',dpi=200)
print('PASS: finite curve data, speed and length-bound figure saved; quadrature uses 2400 intervals.')

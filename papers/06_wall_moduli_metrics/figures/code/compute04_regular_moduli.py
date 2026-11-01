from pathlib import Path
import json
import sympy as S

x,y,M,rho,tau=S.symbols('x y M rho tau',real=True)
# Raw monomials suffice for relation/Jacobian ranks. Orthonormal monomials
# are a diagonal invertible coordinate rescaling, used for norm identities.
basis3=[x**3,x*x*y,x*y*y,y**3]
basis1=[x,y]
variables=S.symbols('a0:8 c0:8 b0:4 d0:4')
a=variables[:8];c=variables[8:16];b=variables[16:20];d=variables[20:24]
poly=lambda coeff,basis:sum(i*j for i,j in zip(coeff,basis))
A=[poly(a[:4],basis3),poly(a[4:],basis3)]
C=[poly(c[:4],basis3),poly(c[4:],basis3)]
B=[poly(b[:2],basis1),poly(b[2:],basis1)]
D=[poly(d[:2],basis1),poly(d[2:],basis1)]
O=S.expand(sum(A[i]*B[i]-C[i]*D[i] for i in range(2)))
coeff=S.Matrix([O.coeff(x,4-i).coeff(y,i) for i in range(5)])
Jac=coeff.jacobian(variables)
point=[1,0,0,0,M,0,0,1]+[1,0,0,0,0,0,0,0]+[1,0,0,0]+[1,0,0,0]
sub=dict(zip(variables,point));assert coeff.subs(sub)==S.zeros(5,1)
J=Jac.subs(sub)
# a0 cubic variations give first four quartics; b1.y supplies y^4.
minor=S.factor(J[:,[0,1,2,3,19]].det());assert minor==1
old=point.copy();old[7]=0
assert Jac.subs(dict(zip(variables,old))).rank()==4
v=S.zeros(24,1);v[5]=1
assert J*v==S.zeros(5,1)
# Action of three zero-sum diagonal generators, written in all four
# coordinates first. Coefficient-space adjoint below uses orthonormal bases.
Inc=S.Matrix([[1,1,0,0],[-1,0,1,0],[0,-1,0,1],[0,0,-1,-1]])
T=S.zeros(24,4)
for edge,(start,stop) in enumerate([(0,8),(8,16),(16,20),(20,24)]):
 for i in range(start,stop):
  for j in range(4):T[i,j]=point[i]*Inc[j,edge]
assert J*T==S.zeros(5,4)
assert T.rank()==3
assert v.T*T==S.zeros(1,4)
assert (S.Matrix(point[:8]).dot(S.Matrix(point[:8])))==M*M+2
# The independently normalized cubic monomials all have unit norm.
# Hence the transverse tangent has norm one and its balancing factor
# squared is w12/(M^2+2), without a finite horizontal correction.
weights=S.diag(1-tau,tau,S.Rational(3,5)-tau,tau)
L=Inc*weights*Inc.T
weak_vector=S.Matrix([-1,-1,3,-1])
assert S.simplify((weak_vector.T*L*weak_vector)[0]/12)==S.Rational(8,3)*tau
balance=S.expand(rho*tau*tau-(1-tau)*(S.Rational(3,5)-tau))
# Parametrize rho by tau to check unrestricted large-rho asymptotics.
rho_of_tau=(1-tau)*(S.Rational(3,5)-tau)/(tau*tau)
assert S.simplify(balance.subs(rho,rho_of_tau))==0
assert S.limit(rho_of_tau*tau*tau,tau,0,dir='+')==S.Rational(3,5)
assert S.limit(1-tau,tau,0)==1
out={'arithmetic':'exact SymPy','relation_vanishes':True,'regular_jacobian_minor':str(minor),'regular_rank':5,'old_example_rank':4,'complex_orbit_rank':3,'quotient_complex_dimension':24-5-3,'transverse_tangent_closed':True,'transverse_tangent_orthogonal_to_complex_orbit':True,'rho':'M^2 + 2','finite_tangent_squared_norm':'(l - delta)/(M^2 + 2)','weak_rayleigh_quotient':'8 delta / 3','rho_tau_squared_limit':'3/5','scope':'Finite algebra only; not the PDE, Hodge projection or gauge-quotient proof.'}
(Path(__file__).resolve().parents[1]/'output'/'regular_moduli_results.json').write_text(json.dumps(out,indent=2)+'\n')
print('All exact relation, regularity, orbit, tangent and scale assertions passed.')

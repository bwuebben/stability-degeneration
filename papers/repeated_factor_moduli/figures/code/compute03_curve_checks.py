"""Checks of the finite curve formulas; no elliptic or metric theorem."""
from pathlib import Path
import json
import numpy as np
from finite_models import curve_data, eigenvalues, MASS, balanced_rows

xi=np.array([1,1,0,0,-4,1])
assert xi@MASS@xi == 20
assert np.linalg.norm(balanced_rows(0)@xi) == 0
assert np.linalg.matrix_rank(balanced_rows(0)) == 5
assert np.linalg.matrix_rank(balanced_rows(.2)) == 6
ratios=[eigenvalues(e)[0]/e**2 for e in [1e-2,1e-3,1e-4]]
assert abs(ratios[-1]-1.25)<1e-6
records=[]
for t in [.01,.1,.2]:
    for eps in [1e-4,.05,.2,.45]:
        f,p,lower,params,x=curve_data(t,eps)
        f0,p0,_,_,_=curve_data(t,eps,rotated=False)
        assert f >= lower*(1-1e-10)
        assert f0 < 1e-18
        assert np.allclose(p,p0)
        records.append({"t":t,"epsilon":eps,"finite_bundle_squared":f,
                        "polarization_squared":p,"orthogonal_lower":lower,
                        "unrotated_bundle_squared":f0})
out={"scope":"finite numerical and exact matrix checks; no HE computation",
     "boundary_mass":20,"ranks":[5,6],"weak_ratios":ratios,
     "curve_checks":records,"pass":True}
path=Path(__file__).resolve().parents[1]/"output"/"curve_checks.json"
path.write_text(json.dumps(out,indent=2)+"\n")
print(json.dumps({"pass":True,"curve_checks":len(records),
                  "weak_ratios":ratios},indent=2))

#!/usr/bin/env python3
"""Independent Riemann--Roch arithmetic for the three Higgs strata.

Uses separately proved stable-factor Hom vanishings and the stated SO(32)
representation decomposition. These integers alone do not verify the complete
physical fluctuation complex or a quantum anomaly calculation.
"""
import json


def h1_bundle(rank, c2):
    return c2 - 2*rank


def neutral_hypers(rank, c2):
    return rank*c2-rank*rank+1


def h1_tensor(r1, k1, r2, k2):
    return r2*k1+r1*k2-2*r1*r2


def h1_exterior(rank, c2):
    return (rank-2)*c2-rank*(rank-1)


def main():
    generic = {
        "E_neutral": neutral_hypers(6, 14),
        "B_neutral": neutral_hypers(3, 10),
        "exterior_E": h1_exterior(6, 14),
        "exterior_B": h1_exterior(3, 10),
        "E_tensor_B": h1_tensor(6, 14, 3, 10),
        "E_tensor_B_dual": h1_tensor(6, 14, 3, 10),
        "SO14_E": 14*h1_bundle(6, 14),
        "SO14_B": 14*h1_bundle(3, 10),
        "geometric": 20,
    }
    axis = dict(generic)
    # Two nonisomorphic rank-three summands: 2*13 neutral hypers and
    # 24 off-diagonal charged hypers in place of the stable E's 49.
    axis["E_neutral"] = 2*neutral_hypers(3, 7)
    axis["E_offdiagonal"] = h1_tensor(3, 7, 3, 7)
    assert axis["exterior_E"] == 2*h1_exterior(3, 7)+h1_tensor(3, 7, 3, 7)
    enhanced = dict(generic)
    enhanced["E_neutral"] = 4*neutral_hypers(3, 7)
    tables = {"enhanced": (enhanced, 96), "axis": (axis, 94), "generic": (generic, 93)}
    result = {}
    for label, (table, n_v) in tables.items():
        n_h = sum(table.values())
        assert n_h-n_v == 244
        assert n_h-n_v+29 == 273
        result[label] = {"hypermultiplet_dimension_table": table, "n_H": n_h, "n_V": n_v}
    assert [result[k]["n_H"] for k in ("enhanced", "axis", "generic")] == [340, 338, 337]
    # Failure control: keeping the enhanced hyper count after Higgsing fails.
    assert result["enhanced"]["n_H"]-result["generic"]["n_V"] != 244
    print(json.dumps({"checks": "passed", "strata": result,
                      "scope": "Riemann--Roch and representation arithmetic, not a full mass-spectrum proof"}, indent=2))


if __name__ == "__main__":
    main()

"""Exact finite checks of the several-factor confinement identities.

This program checks graph algebra and a stated remainder bound. It does
not prove the invariant-manifold theorem or an actual HYM flow theorem."""

from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json


def connected(n, edges):
    seen = {0}
    while True:
        bigger = seen | {
            v for edge in edges if any(w in seen for w in edge) for v in edge
        }
        if bigger == seen:
            return len(seen) == n
        seen = bigger


def supports(n, p):
    all_edges = [(i, j) for i in range(p) for j in range(p, n)]
    for mask in range(1, 1 << len(all_edges)):
        edges = [e for z, e in enumerate(all_edges) if (mask >> z) & 1]
        if len(edges) >= n - 1 and connected(n, edges):
            yield edges


def incidence(n, edges, weights):
    ans = [F(0) for _ in range(n)]
    for (i, j), k in zip(edges, weights):
        ans[i] += k
        ans[j] -= k
    return ans


def inv_mass_square(x, ranks):
    return sum((z * z / r for z, r in zip(x, ranks)), F(0))


def derivative(edges, weights, velocity):
    return 2 * sum(
        (k * (velocity[i] - velocity[j]) for (i, j), k in zip(edges, weights)),
        F(0),
    )


def centered_direction(n):
    raw = [F(i - (n - 1), 1) for i in range(n)]
    avg = sum(raw) / n
    return [z - avg for z in raw]


def main():
    sigma = F(1, 10_000)
    sqrt_sigma = F(1, 100)
    counts = {
        "support_graphs": 0,
        "rank_weight_cases": 0,
        "perturbed_confinement_cases": 0,
    }
    by_n = {}
    smallest = F(1)
    largest_weight_ratio = F(1)
    for n in range(2, 7):
        by_n[str(n)] = 0
        for p in range(1, n):
            for edges in supports(n, p):
                counts["support_graphs"] += 1
                by_n[str(n)] += 1
                ec = len(edges)
                raw_patterns = [
                    [F(1) for _ in edges],
                    [F(1, z + 1) for z in range(ec)],
                    [
                        F(1, 10 ** ((80 * z) // max(1, ec - 1)))
                        for z in range(ec)
                    ],
                ]
                for ranks, raw in product(
                    [
                        [F(1) for _ in range(n)],
                        [F(1 + (v % 3)) for v in range(n)],
                    ],
                    raw_patterns,
                ):
                    k = [sigma * z / sum(raw) for z in raw]
                    smallest = min(smallest, *k)
                    largest_weight_ratio = max(largest_weight_ratio, max(k) / min(k))
                    flow = incidence(n, edges, k)
                    ri = sum(ranks[:p])
                    rj = sum(ranks[p:])
                    rho = 1 / ri + 1 / rj
                    cp2 = 1 / min(ranks[:p]) + 1 / min(ranks[p:])
                    norm_flow2 = inv_mass_square(flow, ranks)
                    assert sum(flow) == 0
                    assert sum(flow[:p]) == sigma
                    assert sum(flow[p:]) == -sigma
                    assert rho * sigma**2 <= norm_flow2 <= cp2 * sigma**2
                    counts["rank_weight_cases"] += 1

                    q = centered_direction(n)
                    qnorm2 = inv_mass_square(q, ranks)
                    small = rho / (16 * (1 + cp2) * (1 + qnorm2))
                    eps = rho / (8 * cp2)
                    forcing_patterns = [
                        [F(0) for _ in flow],
                        [eps * z for z in flow],
                        [-eps * z for z in flow],
                        [small * sigma * z for z in q],
                    ]
                    for forcing, sign in product(forcing_patterns, [-1, 0, 1]):
                        # E is a rank-weighted trace-zero cubic remainder.
                        remainder = [
                            sign * small * sigma * sqrt_sigma * z / r
                            for z, r in zip(q, ranks)
                        ]
                        assert sum(r * e for r, e in zip(ranks, remainder)) == 0
                        beta2 = inv_mass_square(forcing, ranks)
                        assert beta2 <= rho**2 * sigma**2 / (16 * cp2)
                        # |E|_M <= K sigma^(3/2), checked without irrational roots.
                        kbound2 = small**2 * qnorm2
                        enorm2 = sum(
                            (r * e * e for r, e in zip(ranks, remainder)), F(0)
                        )
                        assert enorm2 <= kbound2 * sigma**3
                        assert 4 * cp2 * kbound2 * sigma <= rho**2
                        velocity = [
                            (b - z) / r + e
                            for b, z, r, e in zip(forcing, flow, ranks, remainder)
                        ]
                        exact = derivative(edges, k, velocity)
                        cross = sum(
                            (z * b / r for z, b, r in zip(flow, forcing, ranks)),
                            F(0),
                        )
                        error = sum(
                            (z * e for z, e in zip(flow, remainder)), F(0)
                        )
                        assert exact == -2 * norm_flow2 + 2 * cross + 2 * error
                        assert abs(2 * cross) <= rho * sigma**2 / 2
                        assert abs(2 * error) <= rho * sigma**2
                        assert exact <= -rho * sigma**2 / 2
                        counts["perturbed_confinement_cases"] += 1

    # A weak edge can grow while the total conductance still decreases.
    eps = F(1, 10**6)
    edges = [(0, 2), (0, 3), (1, 2), (1, 3)]
    k = [eps, eps**2, eps**2, F(1)]
    flow = incidence(4, edges, k)
    forcing = [2 * eps, eps, -2 * eps, -eps]
    velocity = [b - z for b, z in zip(forcing, flow)]
    edge_derivatives = [
        2 * w * (velocity[i] - velocity[j])
        for (i, j), w in zip(edges, k)
    ]
    assert edge_derivatives[0] > 0
    assert edge_derivatives[-1] < 0
    assert sum(edge_derivatives) < 0
    assert inv_mass_square(forcing, [F(1)] * 4) <= sum(k)**2 / 32

    # Controls: reversing the dynamical sign, a missing factor two, an
    # oversized slope demand, and a non-bipartite directed cycle.
    e = [(0, 1)]
    w = [F(1, 100)]
    v = incidence(2, e, w)
    good = [-z for z in v]
    assert derivative(e, w, good) == -2 * inv_mass_square(v, [F(1), F(1)])
    assert derivative(e, w, good) != -inv_mass_square(v, [F(1), F(1)])
    assert derivative(e, w, v) > 0
    oversized_b = [2 * z for z in v]
    assert derivative(e, w, [b - z for b, z in zip(oversized_b, v)]) > 0
    cyc = [(1, 0), (2, 1), (0, 2)]
    cyc_flow = incidence(3, cyc, [F(1)] * 3)
    assert cyc_flow == [0, 0, 0]
    assert inv_mass_square(cyc_flow, [F(1)] * 3) < F(9)

    result = {
        "scope": "Exact finite graph algebra with a stated cubic remainder; no invariant-manifold or PDE verification",
        "support_graphs_by_vertex_count": by_n,
        **counts,
        "weight_patterns_per_support": 3,
        "rank_patterns_per_support": 2,
        "largest_exact_weight_ratio": str(largest_weight_ratio),
        "smallest_weight_approximation": f"{float(smallest):.6e}",
        "mixed_edge_velocity_with_decreasing_total_detected": True,
        "wrong_dynamical_sign_detected": True,
        "missing_factor_two_detected": True,
        "oversized_parameter_demand_detected": True,
        "non_bipartite_coercivity_failure_detected": True,
    }
    output = (Path(__file__).resolve().parents[1] / "data/graph_dynamics_checks.json")
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

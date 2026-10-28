"""Exact finite checks of the order-ideal, Harder–Narasimhan and four-cycle formulas.

Uses only the Python standard library. These computations do not prove
the geometric or infinite-dimensional assertions of the paper."""

from fractions import Fraction as F
from itertools import combinations, permutations, product
from collections import Counter
import json
from pathlib import Path


VERTICES = frozenset(range(4))
EDGES = ((0, 2), (0, 3), (1, 2), (1, 3))
IDEALS = tuple(
    frozenset(subset)
    for size in range(5)
    for subset in combinations(VERTICES, size)
    if all(q not in subset or p in subset for p, q in EDGES)
)


def slope(z, lower, upper):
    return sum(z[p] for p in upper - lower) / F(len(upper - lower))


def greedy(z):
    result = [frozenset()]
    while result[-1] != VERTICES:
        previous = result[-1]
        candidates = [j for j in IDEALS if previous < j]
        maximum = max(slope(z, previous, j) for j in candidates)
        maximizers = [j for j in candidates if slope(z, previous, j) == maximum]
        largest = frozenset().union(*maximizers)
        assert largest in maximizers
        result.append(largest)
    return tuple(result)


def flags(previous=frozenset()):
    if previous == VERTICES:
        yield (previous,)
        return
    for j in IDEALS:
        if previous < j:
            for continuation in flags(j):
                yield (previous,) + continuation


ALL_FLAGS = tuple(flags())


def hn_characterization(z, flag):
    quotients = [slope(z, a, b) for a, b in zip(flag, flag[1:])]
    if any(a <= b for a, b in zip(quotients, quotients[1:])):
        return False
    return all(
        slope(z, a, k) <= slope(z, a, b)
        for a, b in zip(flag, flag[1:])
        for k in IDEALS
        if a < k < b
    )


def poly_add(a, b):
    c = [F(0)] * max(len(a), len(b))
    for i, x in enumerate(a):
        c[i] += x
    for i, x in enumerate(b):
        c[i] += x
    return c


def poly_mul(a, b):
    c = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i + j] += x * y
    return c


def characteristic_polynomial(matrix):
    result = [F(0)]
    for perm in permutations(range(4)):
        inversions = sum(perm[i] > perm[j] for i in range(4) for j in range(i + 1, 4))
        term = [F((-1) ** inversions)]
        for i, j in enumerate(perm):
            term = poly_mul(term, [-matrix[i][j], F(i == j)])
        result = poly_add(result, term)
    return result


def laplacian(a, b):
    weights = (a, b, b, a)
    result = [[F(0)] * 4 for _ in range(4)]
    for (p, q), weight in zip(EDGES, weights):
        result[p][p] += weight
        result[q][q] += weight
        result[p][q] -= weight
        result[q][p] -= weight
    return result


def main():
    checked = 0
    observed = set()
    witnesses = {}
    for first in product(range(-2, 3), repeat=3):
        z = tuple(map(F, (*first, -sum(first))))
        expected = [flag for flag in ALL_FLAGS if hn_characterization(z, flag)]
        assert len(expected) == 1
        obtained = greedy(z)
        assert expected == [obtained]
        observed.add(obtained)
        witnesses.setdefault(obtained, z)
        checked += 1

    # Mandatory orientation control: upward cone has negative minimal demands.
    interior = tuple(map(F, (-1, -1, 1, 1)))
    assert greedy(interior) == (frozenset(), VERTICES)
    assert greedy(tuple(-x for x in interior)) == (
        frozenset(), frozenset((0, 1)), VERTICES
    )

    spectra = []
    for q in (F(1, 3), F(1, 2), F(1), F(2), F(3), F(10)):
        a, b = q / (1 + q), 1 / (1 + q)
        eigenvalues = [F(0), 2 * a, 2 * b, F(2)]
        prediction = [F(1)]
        for eigenvalue in eigenvalues:
            prediction = poly_mul(prediction, [-eigenvalue, F(1)])
        assert characteristic_polynomial(laplacian(a, b)) == prediction
        assert a * a / (b * b) == q * q
        spectra.append({"sqrt_cross_ratio": str(q), "eigenvalues": list(map(str, eigenvalues))})

    # Failure control: vanishing alternating edges disconnects the graph.
    disconnected = characteristic_polynomial(laplacian(F(0), F(1)))
    assert disconnected[0] == disconnected[1] == 0
    assert disconnected[2] != 0

    def components(subset):
        remaining = set(subset)
        answer = 0
        while remaining:
            pending = [remaining.pop()]
            answer += 1
            while pending:
                vertex = pending.pop()
                adjacent = {
                    q if p == vertex else p
                    for p, q in EDGES if p == vertex or q == vertex
                } & remaining
                remaining -= adjacent
                pending.extend(adjacent)
        return answer

    dimensions = Counter()
    witness_records = []
    for flag in ALL_FLAGS:
        dimension = 3 - sum(components(b - a) - 1 for a, b in zip(flag, flag[1:]))
        dimensions[dimension] += 1
        witness_records.append({
            "flag": [sorted(i) for i in flag],
            "dimension_in_Z": dimension,
            "exact_demand_witness": list(map(str, witnesses[flag])),
        })
    assert dimensions == {1: 1, 2: 4, 3: 13}

    result = {
        "status": "exact finite verification; not a geometric or analytic proof",
        "crown_ideals": [sorted(i) for i in IDEALS],
        "proper_nonempty_ideal_count": len(IDEALS) - 2,
        "all_ideal_flags": len(ALL_FLAGS),
        "integer_zero_sum_demands_checked": checked,
        "distinct_hn_flags_observed": len(observed),
        "hn_flag_dimensions": dict(sorted(dimensions.items())),
        "hn_witnesses": witness_records,
        "greedy_vs_exhaustive_flags": "PASS",
        "orientation_failure_control": "PASS",
        "spectrum_determinant_comparisons": spectra,
        "disconnected_failure_control": "PASS; two zero modes",
        "square_slice_vertices": [
            [-1, 0, 1, 0], [-1, 0, 0, 1], [0, -1, 1, 0], [0, -1, 0, 1]
        ],
    }
    target = Path(__file__).resolve().parents[1] / "data/topology_constraints_checks.json"
    target.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k not in ("spectrum_determinant_comparisons", "hn_witnesses")}, indent=2))


if __name__ == "__main__":
    main()

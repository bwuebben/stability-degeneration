"""Exact graph-lattice calculations for the integral-monodromy example.

These calculations do not prove geometric realization or analytic recovery.
Run with: uv run --with sympy python <this-file>
Adapted from the accepted R11 computation without changing its graph data.
"""
from pathlib import Path
from collections import deque
from itertools import combinations
from math import gcd
import json
import sympy as s
from sympy.matrices.normalforms import smith_normal_form


def graph_record(name, edges):
    n = 6
    edges = sorted(tuple(sorted(e)) for e in edges)
    assert len(set(edges)) == len(edges) == 9
    incidence = s.zeros(n, len(edges))
    tree, chords = [], []
    forest = [{i} for i in range(n)]
    for k, (i, j) in enumerate(edges):
        incidence[i, k], incidence[j, k] = -1, 1
        ci = next(c for c in forest if i in c)
        cj = next(c for c in forest if j in c)
        if ci is cj:
            chords.append(k)
        else:
            tree.append(k)
            forest.remove(ci)
            forest.remove(cj)
            forest.append(ci | cj)
    assert len(forest) == 1 and len(tree) == 5 and len(chords) == 4
    assert all(sum(abs(incidence[i, j]) for j in range(9)) == 3 for i in range(n))
    adjacency = [[] for _ in range(n)]
    for k in tree:
        i, j = edges[k]
        adjacency[i].append((j, k, 1))
        adjacency[j].append((i, k, -1))
    cycles = s.zeros(9, 4)
    for col, chord in enumerate(chords):
        i, j = edges[chord]
        cycles[chord, col] = 1
        queue = deque([(j, [])])
        visited = {j}
        while queue:
            vertex, path = queue.popleft()
            if vertex == i:
                for k, sign in path:
                    cycles[k, col] = sign
                break
            for nxt, k, sign in adjacency[vertex]:
                if nxt not in visited:
                    visited.add(nxt)
                    queue.append((nxt, path + [(k, sign)]))
        else:
            raise AssertionError('Missing tree path')
    assert incidence * cycles == s.zeros(6, 4)
    assert cycles.extract(chords, list(range(4))) == s.eye(4)
    gram = cycles.T * cycles
    laplacian = incidence * incidence.T
    tree_count = int(laplacian[:5, :5].det())
    assert gram.det() == tree_count
    smith = smith_normal_form(gram, domain=s.ZZ)
    diagonal = [abs(int(smith[i, i])) for i in range(4)]
    assert s.prod(diagonal) == tree_count
    minor_gcds = []
    for size in range(1, 5):
        divisor = 0
        for rows in combinations(range(4), size):
            for cols in combinations(range(4), size):
                divisor = gcd(divisor, int(gram.extract(rows, cols).det()))
        minor_gcds.append(divisor)
    assert minor_gcds == [int(s.prod(diagonal[:k])) for k in range(1, 5)]
    nilpotent = s.zeros(8)
    nilpotent[:4, 4:] = gram
    symplectic = s.zeros(8)
    symplectic[:4, 4:] = s.eye(4)
    symplectic[4:, :4] = -s.eye(4)
    monodromy = s.eye(8) + nilpotent
    assert nilpotent**2 == s.zeros(8)
    assert monodromy.T * symplectic * monodromy == symplectic
    return {'graph': name, 'edges': edges, 'tree_edge_indices': tree,
            'cycle_basis': [list(map(int, cycles.row(i))) for i in range(9)],
            'cycle_gram': [list(map(int, gram.row(i))) for i in range(4)],
            'smith_diagonal': diagonal, 'minor_gcds': minor_gcds,
            'spanning_trees': tree_count,
            'genus_if_all_components_rational': 4,
            'nilpotent_rank': nilpotent.rank()}


def main():
    prism = [(0, 1), (1, 2), (0, 2), (3, 4), (4, 5), (3, 5),
             (0, 3), (1, 4), (2, 5)]
    bipartite = [(i, j) for i in range(3) for j in range(3, 6)]
    records = [graph_record('triangular_prism', prism),
               graph_record('K3_3', bipartite)]
    assert records[0]['smith_diagonal'] == [1, 1, 5, 15]
    assert records[1]['smith_diagonal'] == [1, 3, 3, 9]
    assert records[0]['minor_gcds'] == [1, 1, 5, 75]
    assert records[1]['minor_gcds'] == [1, 3, 9, 81]
    assert records[0]['cycle_gram'] == [
        [3, -1, 1, 2], [-1, 4, 1, -2], [1, 1, 4, 2], [2, -2, 2, 5]]
    assert records[1]['cycle_gram'] == [
        [4, 2, 2, 1], [2, 4, 1, 2], [2, 1, 4, 2], [1, 2, 2, 4]]
    result = {'scope': 'Exact finite integral lattice calculations only; geometric realization and recovery from analytic observations require separate proofs.',
              'records': records}
    output = Path(__file__).resolve().parents[1] / 'output' / 'integral_graph_pairs.json'
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + '\n')
    for r in records:
        print(r['graph'], 'minor gcds:', r['minor_gcds'],
              'Smith:', r['smith_diagonal'], 'det:', r['spanning_trees'])
    print('PASS: cycle lattices, minor gcds, Smith forms and symplectic transvections.')


if __name__ == '__main__':
    main()

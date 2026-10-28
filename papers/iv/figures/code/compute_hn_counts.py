#!/usr/bin/env python3
"""Compare direct ideal-chain counts with ordered-partition predictions.

Only exact stdlib arithmetic is used. The graph is K_(2,2k), k=1,2,3.
This finite check does not prove the all-k geometric/HN statements.
"""
from pathlib import Path
from math import comb
import json


def ordered_bell(n):
    values = [1]
    for m in range(1, n + 1):
        values.append(sum(comb(m, j) * values[m-j] for j in range(1, m+1)))
    return values[-1]


def ideal_chain_count(a, b):
    total = (1 << (a+b)) - 1
    minima = (1 << a) - 1
    ideals = [mask for mask in range(total+1)
              if not (mask >> a) or (mask & minima) == minima]
    ways = {0: 1}
    for mask in sorted(ideals[1:], key=int.bit_count):
        ways[mask] = sum(ways[sub] for sub in ways if sub != mask and sub & mask == sub)
    return len(ideals), ways[total]


def main():
    records = []
    for k, expected in ((1, 18), (2, 450), (3, 28098)):
        ideals, count = ideal_chain_count(2, 2*k)
        prediction = 2 * ordered_bell(2) * ordered_bell(2*k)
        assert count == prediction == expected
        assert ideals == 2**(2*k) + 3
        records.append({'k': k, 'ideal_count': ideals, 'HN_flags': count,
                        'ordered_Bell_2k': ordered_bell(2*k)})
    result = {'status': 'PASS; exact finite check', 'records': records,
              'method': 'dynamic count of ideal chains versus ordered partitions'}
    path = Path(__file__).resolve().parents[1] / 'data/hn_counts.json'
    path.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()

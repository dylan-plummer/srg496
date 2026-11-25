"""Intrinsic fibres: 4-sets of vertices whose adjacency rows sum (mod 2) to 0 or to the
all-ones vector, outside the 4-set itself.  Found by hashing the XOR of pairs of rows of
A + I (a 4-set {a,b,c,d} qualifies iff row_a^row_b == row_c^row_d or its complement,
up to the four coordinates in the set).  Reports how many such 4-sets lie inside one fibre
(16 special columns, then 60 fibres of 8 columns)."""
import sys
from collections import defaultdict
from itertools import combinations
from pathlib import Path
import numpy as np

def fibre(v):
    return -1 if v < 16 else (v - 16) // 8

def zero_sum_quads(A):
    n = len(A)
    M = (A + np.eye(n, dtype=np.int64)) % 2
    rows = [int(''.join(map(str, r)), 2) for r in M]
    full = (1 << n) - 1
    buckets = defaultdict(list)
    for a, b in combinations(range(n), 2):
        x = rows[a] ^ rows[b]
        key = min(x, x ^ full)
        buckets[key].append((a, b))
    quads = set()
    for key, prs in buckets.items():
        for (a, b), (c, d) in combinations(prs, 2):
            s = {a, b, c, d}
            if len(s) == 4:
                quads.add(tuple(sorted(s)))
    return quads

def main(path):
    from independent_rebuild import read_g6
    A = read_g6(Path(path))
    Q = zero_sum_quads(A)
    inside = [q for q in Q if len({fibre(v) for v in q}) == 1]
    print(f'{path}: {len(Q)} zero-sum 4-sets (exact, up to complement), {len(inside)} inside one fibre')
    from collections import Counter
    print('fibres hit:', Counter(fibre(q[0]) for q in inside).most_common(3), '...',
          'distinct fibres', len({fibre(q[0]) for q in inside}))
    print('non-fibre examples:', [q for q in Q if q not in set(inside)][:5])

if __name__ == '__main__':
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'audit'))
    main(sys.argv[1] if len(sys.argv) > 1 else 'graphs/srg_496_231_102_112.g6')

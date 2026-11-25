"""The 496 frame is half of a Steiner equiangular tight frame of the S(2,5,61).

Steiner ETF (Fickus-Mixon-Tremain): rows = the 183 blocks, and at each point p the 16 columns
y in F_2^4 have entry (-1)^{lambda_p(B).y} on block B through p, where lambda_p is a bijection
from the 15 blocks through p to the nonzero characters of F_2^4.
Labeling at an ordinary point p in part P: lambda_p(P) = (0,0,0,1); in a tile T where p is
the K6 edge ab with character c and root sign s, star a gets (c,0) if s=-1 else (c,1), and
star b the other one.  Then the eight columns with y_4 = 1, negated, are exactly the eight
columns of p in the 496 frame.  At the added point, lambda(part j) = j+1 gives the 16 special
columns.  The script checks both statements and that the full 976 columns form an ETF of
rank 183 (Seidel eigenvalues 65, -15) -- the Steiner two-graph covered by the corollary in the
442 paper -- while the 496-column half has rank 155.
"""
import os
os.environ.setdefault('OMP_NUM_THREADS', '8')
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from independent_rebuild import (tiles_from, read_starter, described_labels, frame, part, chi)

def main():
    tiles = tiles_from(read_starter())
    chars, signs = described_labels(tiles)
    F = frame(tiles, chars, signs)
    nrows = F.shape[0]
    cols = []
    # added point: 16 columns
    for y in range(16):
        v = np.zeros(nrows, dtype=np.int64)
        for j in range(15):
            v[j] = chi(j + 1, y)
        cols.append(v)
    halves = []
    for p in range(60):
        mine = [t for t in range(28) if p in tiles[t]]
        lab = {('part',): 0b1000}            # character (0,0,0,1) as bit 3
        for t in mine:
            a, b = tiles[t][p]
            c, s = chars[p][t], signs[p][t]
            lab[(t, a)] = c | (0 if s == -1 else 0b1000)
            lab[(t, b)] = c | (0b1000 if s == -1 else 0)
        assert sorted(lab.values()) == list(range(1, 16))
        for y in range(16):
            v = np.zeros(nrows, dtype=np.int64)
            v[part(p)] = chi(lab[('part',)], y)
            for t in mine:
                a, b = tiles[t][p]
                v[15 + 6 * t + a] = chi(lab[(t, a)], y)
                v[15 + 6 * t + b] = chi(lab[(t, b)], y)
            cols.append(v)
            if y & 0b1000:
                halves.append(-v)
    Phi = np.array(cols).T
    sub = np.concatenate([F[:, :16], np.array(halves).T], axis=1)
    print('496 frame = special columns + negated y4=1 halves of the Steiner columns:',
          np.array_equal(sub, F))
    G = Phi.T @ Phi
    S = G - 15 * np.eye(976, dtype=np.int64)
    off = S[~np.eye(976, dtype=bool)]
    print('Steiner frame: 976 columns, equiangular:', bool(np.all(np.abs(off) == 1)),
          ' tight:', np.array_equal(Phi @ Phi.T, 80 * np.eye(nrows, dtype=np.int64)))
    ev = np.round(np.linalg.eigvalsh(S.astype(float))).astype(int)
    v, m = np.unique(ev, return_counts=True)
    print('Steiner Seidel spectrum', dict(zip(v.tolist(), m.tolist())))
    print('rank of the 496 half:', np.linalg.matrix_rank(F.astype(float)))

if __name__ == '__main__':
    main()

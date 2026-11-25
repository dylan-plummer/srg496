"""Independent check of the Hadamard half-frames and the 976-point completion (written from the text only).

  * The five normalized Hadamard matrices of order 16 in hadamards.json (row r, column x entry
    1 - 2*bit_x(word_r), first row all ones): H H^T = 16 I; on the 16 columns, the four-sets whose
    product over the 15 other rows is constant number 140, 76, 44, 28, 28; the group of column
    permutations preserving the 15 unordered pairs of half-blocks has orders 322560, 9216, 1536, 2688, 2688,
    and the five structures are pairwise non-isomorphic (nauty).
  * Replacing the Walsh columns at infinity by any of them: random admissible members of all three designs give
    srg(496,231,102,112), srg(496,240,120,112), srg(495,238,109,119), and after a random switching and
    relabeling the zero-sum four-sets of rows of S - I number 840 + (140, 76, 44, 28, 28).
  * The 976-point completion: a second group of eight columns at each point, -1 on the part row and
    chi_d(x)(e_a + e_b) in each tile, with an independent random labeling. Then F F^T = 80 I, the Seidel
    matrix has eigenvalues 65 and -15, the cell-constant kernel word gives srg(976,480,248,224), whose
    first 496 vertices induce the srg(496,240,120,112) of the half-frame, and a descendant is
    srg(975,462,205,231).
  * The lower bounds: 486 times the earlier bound, and ceil(30^60 / 84).
"""
import json
import random
import sys
from fractions import Fraction
from itertools import combinations
from math import ceil, factorial
from pathlib import Path
import numpy as np
import pynauty

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parent.parent / 'audit'))
sys.path.insert(0, str(HERE.parent.parent / 'family'))
from designs_audit import random_decoration, check_fibres, STARTERS
from independent_rebuild import tiles_from, read_starter, frame, regular_graph, described_z2, check_srg, chi, part


def matrices():
    recs = json.loads((HERE / 'hadamards.json').read_text())
    Hs = []
    for r in recs:
        H = np.array([[1 - 2 * ((w >> x) & 1) for x in range(16)] for w in r['matrix']], dtype=np.int64)
        assert np.all(H[0] == 1) and np.array_equal(H @ H.T, 16 * np.eye(16, dtype=np.int64))
        Hs.append(H)
    return Hs


def special_planes(H):
    n = 0
    for q in combinations(range(16), 4):
        prod = H[1:, q[0]] * H[1:, q[1]] * H[1:, q[2]] * H[1:, q[3]]
        n += bool(np.all(prod == prod[0]))
    return n


def pair_structure(H):
    """Columns 0..15, half-blocks 16..45, row pairs 46..60."""
    adj = {v: set() for v in range(61)}
    for r in range(1, 16):
        for k, sgn in enumerate((1, -1)):
            h = 16 + 2 * (r - 1) + k
            adj[h] |= {x for x in range(16) if H[r, x] == sgn} | {46 + r - 1}
    for v in list(adj):
        for w in adj[v]:
            adj[w].add(v)
    g = pynauty.Graph(61, adjacency_dict={v: sorted(w) for v, w in adj.items()},
                      vertex_coloring=[set(range(16)), set(range(16, 46)), set(range(46, 61))])
    _, size, exp, _, _ = pynauty.autgrp(g)
    return pynauty.certificate(g), round(size * 10 ** exp)


def install(F, H, rng):
    """Put the 15 non-constant rows of H on the part rows of the 16 columns at infinity, random part order."""
    F = F.copy()
    rows = list(range(1, 16))
    rng.shuffle(rows)
    cols = list(range(16))
    rng.shuffle(cols)
    for g in range(15):
        F[g, :16] = H[rows[g]][cols]
    return F


def second_half(tiles, rng):
    """183 x 480 block: second group of 8 columns at each point."""
    G = np.zeros((183, 480), dtype=np.int64)
    for p in range(60):
        mine = [t for t in range(28) if p in tiles[t]]
        cs = list(range(1, 8))
        rng.shuffle(cs)
        for x in range(8):
            col = 8 * p + x
            G[part(p), col] = -1
            for t, c in zip(mine, cs):
                a, b = tiles[t][p]
                G[15 + 6 * t + a, col] += chi(c, x)
                G[15 + 6 * t + b, col] += chi(c, x)
    return G


def kernel_word_976(rng):
    w = np.ones(976, dtype=np.int64)
    for g in range(15):
        cells = [4 * g + i for i in range(4)]
        rng.shuffle(cells)
        for q in cells[:2]:
            w[16 + 8 * q:24 + 8 * q] = -1             # first group
            w[496 + 8 * q:504 + 8 * q] = -1           # second group
    return w


def descendant(S, r):
    T = S * S[r][:, None] * S[r][None, :]
    keep = [v for v in range(len(S)) if v != r]
    return (T[np.ix_(keep, keep)] < 0).astype(np.int64)


def main():
    rng = random.Random(976)
    Hs = matrices()
    planes = [special_planes(H) for H in Hs]
    certs, orders = zip(*(pair_structure(H) for H in Hs))
    print('special four-sets', planes, 'group orders', list(orders))
    assert planes == [140, 76, 44, 28, 28] and list(orders) == [322560, 9216, 1536, 2688, 2688]
    assert len(set(certs)) == 5
    designs = [tiles_from(read_starter(STARTERS / f'design0{i}.starter')) for i in range(3)]
    for k, H in enumerate(Hs):
        for i, tiles in enumerate(designs):
            chars, signs = random_decoration(tiles, rng)
            F = install(frame(tiles, chars, signs), H, rng)
            A, w, S = regular_graph(F, described_z2())
            check_srg(A)
            if i == 0:
                nq, inside = check_fibres(S, rng)
                assert nq == inside == 840 + planes[k], (nq, inside)
            # 976 completion
            Ffull = np.hstack([F, second_half(tiles, rng)])
            gram = Ffull.T @ Ffull
            assert np.array_equal(Ffull @ Ffull.T, 80 * np.eye(183, dtype=np.int64))
            assert np.all(np.diag(gram) == 15) and np.all(np.abs(gram[~np.eye(976, dtype=bool)]) == 1)
            S976 = gram - 15 * np.eye(976, dtype=np.int64)
            wk = kernel_word_976(rng)
            assert not np.any(Ffull @ wk)
            G = ((wk[:, None] * S976 * wk[None, :]) > 0).astype(np.int64)
            check_srg(G, 480, 248, 224)
            G496 = ((wk[:496, None] * S[:, :] * wk[None, :496]) > 0).astype(np.int64)
            check_srg(G496, 240, 120, 112)
            assert np.array_equal(G[:496, :496], G496)
            r = rng.randrange(976)
            check_srg(descendant(S976, r), 462, 205, 231)
            if k == 4 and i == 0:
                r = rng.randrange(496)
                check_srg(descendant(S, r), 238, 109, 119)
        print(f'Hadamard class {k}: 231 rows ok on all three designs; recovery {840 + planes[k]} four-sets; '
              f'976 completion: ETF, srg(976,480,248,224) containing the 496 graph, srg(975,462,205,231)', flush=True)
    old = Fraction(6 ** 56 * 30 ** 4 * factorial(15), 84 * 20160) * 2 ** 152 * (2 ** 22 + 2 * 2 ** 19)
    weights = [Fraction(322560, a) for a in orders]
    assert weights == [1, 35, 210, 120, 120]
    new = old * sum(weights[1:])
    print(f'old bound {float(old):.5e}; new classes {float(new):.5e}; all five {float(old * 486):.5e}; '
          f'parents of one fixed half >= {ceil(Fraction(30 ** 60, 84)):.5e}')


if __name__ == '__main__':
    main()

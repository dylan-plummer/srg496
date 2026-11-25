"""Independent audit of the nine order-84 seeds, the exported graphs and the Burnside count.

  * For each seed make(i, j): the raw Seidel matrix S; nauty computes Aut of the raw graph
    (entries -1), which must have order 84.
  * |Aut(two-graph)| = 84: for one root in every orbit of that group, the descendant's
    nauty certificate and automorphism order; orbit size x stabilizer = 84 and certificates
    differ between orbits.  Across the nine seeds all 126 certificates differ, so the nine
    two-graphs are pairwise non-isomorphic.
  * The exported graphs: srg parameters on all pairs; 56 degree-231 and 126 descendants
    pairwise non-isomorphic, with the stated automorphism orders; every exported 231/240 graph
    lies in one of the nine two-graphs (its descendant at vertex 0 has one of the 126
    certificates) and the 126 descendants are exactly those certificates.
  * Burnside count of the cell-constant kernel words, computed from the nauty group with
    our own fixed-point rule (including the sign of the induced map on characters).
"""
import gzip
import sys
from collections import Counter
from itertools import combinations, product
from math import prod
from pathlib import Path
import numpy as np
import pynauty

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import family as seeds          # the companion constructor of the nine seeds

SRC = HERE / 'graphs'


def g6_decode(line):
    data = line.strip()
    if data[0] == 126:
        n = ((data[1] - 63) << 12) | ((data[2] - 63) << 6) | (data[3] - 63)
        body = data[4:]
    else:
        n, body = data[0] - 63, data[1:]
    bits = np.unpackbits(np.frombuffer(bytes(c - 63 for c in body), dtype=np.uint8)[:, None], axis=1)[:, 2:].ravel()
    A = np.zeros((n, n), dtype=np.int64)
    iu = [(i, j) for j in range(n) for i in range(j)]
    idx = np.array(iu)
    A[idx[:, 0], idx[:, 1]] = bits[:len(iu)]
    return A + A.T


def srg_ok(A, k, lam, mu):
    n = len(A)
    assert np.array_equal(A, A.T) and not A.diagonal().any() and np.all(A.sum(1) == k)
    assert np.array_equal(A @ A, (k - mu) * np.eye(n, dtype=np.int64) + (lam - mu) * A + mu)


def nauty(A):
    g = pynauty.Graph(len(A), adjacency_dict={i: np.nonzero(A[i])[0].tolist() for i in range(len(A))})
    gens, size, exp, orbits, _ = pynauty.autgrp(g)
    return pynauty.certificate(g), round(size * 10 ** exp), gens, orbits


def descendant(S, r):
    T = S * S[r][:, None] * S[r][None, :]
    keep = [v for v in range(len(S)) if v != r]
    return (T[np.ix_(keep, keep)] < 0).astype(np.int64)


def closure(gens, n):
    gens = [tuple(g) for g in gens]
    seen = {tuple(range(n))}
    queue = list(seen)
    for p in queue:
        for g in gens:
            q = tuple(g[i] for i in p)
            if q not in seen:
                seen.add(q)
                queue.append(q)
    return seen


def burnside(F, group):
    """Fixed normalized cell-constant kernel words for each g; returns the two totals x 84."""
    cell = lambda v: (v - 16) // 8                                   # ordinary cells 0..59
    tot = [0, 0]
    for g in group:
        cells = [cell(g[16 + 8 * q]) for q in range(60)]
        assert all(cell(g[16 + 8 * q + x]) == cells[q] for q in range(60) for x in range(8))
        assert all(g[v] < 16 for v in range(16))
        pmap = [cells[4 * h] // 4 for h in range(15)]
        assert all(cells[4 * h + i] // 4 == pmap[h] for h in range(15) for i in range(4))
        # induced map on the special characters, with sign
        sign = {}
        for h in range(15):
            row = F[h, :16]
            moved = np.empty(16, dtype=np.int64)
            moved[[g[x] for x in range(16)]] = row                   # (g.chi)(g x) = chi(x)
            img = pmap[h]
            s = moved * F[img, :16]
            assert abs(s.sum()) == 16, 'special cell is not mapped by an affine map'
            sign[h] = s[0]
        # cycles of parts
        seen, cyc = set(), []
        for h in range(15):
            if h in seen:
                continue
            orbit, x = [], h
            while x not in seen:
                seen.add(x)
                orbit.append(x)
                x = pmap[x]
            L = len(orbit)
            # g^L restricted to the 4 cells of part h
            def gL(q):
                for _ in range(L):
                    q = cells[q]
                return q
            loc = [gL(4 * h + i) - 4 * h for i in range(4)]
            two = sum(1 for a, b in combinations(range(4), 2) if {loc[a], loc[b]} == {a, b})
            fixed = sum(1 for a in range(4) if loc[a] == a)
            cyc.append((orbit, two, fixed))
        star = prod(t for _, t, _ in cyc)
        char = 0
        for i, (orbit, _, fixed) in enumerate(cyc):
            if len(orbit) == 1 and sign[orbit[0]] == 1:
                char += fixed * prod(t for j, (_, t, _) in enumerate(cyc) if j != i)
        tot[0] += star
        tot[1] += char
    return tot


def main():
    all_desc = {}
    totals = []
    for i, j in product(range(3), repeat=2):
        F = seeds.make(i, j)
        S = F.T @ F - 15 * np.eye(496, dtype=np.int64)
        assert np.array_equal((S + 15 * np.eye(496, dtype=np.int64)) @ (S + 15 * np.eye(496, dtype=np.int64)),
                              48 * (S + 15 * np.eye(496, dtype=np.int64)))
        A, _, _ = seeds.switch(F)
        Ex = g6_decode((SRC / f'496_231_{i}{j}.g6').read_bytes())
        _, raw_order, gens, orbits = nauty((S < 0).astype(np.int64) - 0)
        group = closure(gens, 496)
        assert raw_order == len(group) == 84
        assert all(np.array_equal(S[np.ix_(g, g)], S) for g in group)
        reps = sorted({o for o in orbits})
        size = Counter(orbits)
        certs = {}
        for r in reps:
            c, order, _, _ = nauty(descendant(S, r))
            assert size[r] * order == 84
            certs[c] = (r, size[r], order)
        assert len(certs) == len(reps)
        for c in certs:
            assert c not in all_desc
            all_desc[c] = (i, j)
        tot = burnside(F, group)
        assert tot[0] % 84 == 0 and tot[1] % 84 == 0
        totals.append((tot[0] // 84, tot[1] // 84))
        print(f'seed {i}{j}: raw group 84, exported 231 graph equal to switch(make): {np.array_equal(A, Ex)}, '
              f'root orbits {len(reps)} sizes {sorted(size.values())}, two-graph group 84, '
              f'Burnside (constant, character) = {totals[-1]}')
    assert len(all_desc) == 126
    n240 = sum(sum(t) for t in totals)
    print(f'126 distinct descendants over the nine seeds; Burnside total {n240:,} '
          f'(per seed {sum(totals[0]):,})')

    # exported streams
    G231 = [g6_decode(l) for l in gzip.open(SRC / 'all_496_231.g6.gz') if l.strip()]
    G495 = [g6_decode(l) for l in gzip.open(SRC / 'all_495_238.g6.gz') if l.strip()]
    c231 = Counter()
    certs231 = set()
    for A in G231:
        srg_ok(A, 231, 102, 112)
        c, order, _, _ = nauty(A)
        certs231.add(c)
        c231[order] += 1
        S = 1 - np.eye(496, dtype=np.int64) - 2 * A
        assert nauty(descendant(S, 0))[0] in all_desc
    assert len(certs231) == len(G231) == 56
    c495, certs495 = Counter(), set()
    for D in G495:
        srg_ok(D, 238, 109, 119)
        c, order, _, _ = nauty(D)
        certs495.add(c)
        c495[order] += 1
    assert certs495 == set(all_desc)
    print(f'all_496_231: {len(G231)} srg, pairwise non-isomorphic, |Aut| distribution {sorted(c231.items())}, '
          f'all in the nine two-graphs')
    print(f'all_495_238: {len(G495)} srg, = the 126 descendant classes, |Aut| distribution {sorted(c495.items())}')
    for f in sorted(SRC.glob('496_240_*.g6')):
        A = g6_decode(f.read_bytes())
        srg_ok(A, 240, 120, 112)
        S = 2 * A - 1 + np.eye(496, dtype=np.int64)          # the 240 graph is the positive entries
        assert nauty(descendant(S, 0))[0] in all_desc
    for f in sorted(SRC.glob('495_238_*.g6')):
        assert nauty(g6_decode(f.read_bytes()))[0] in all_desc
    print('individual 496_240_*, 495_238_* files: srg, and in the nine two-graphs')


if __name__ == '__main__':
    main()

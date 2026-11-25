"""Independent audit of the three T(6)-tiled S(2,5,61) designs and the lower bound.

Uses only the 15-entry starters and the frame description (audit/independent_rebuild.py).
For each starter:
  * the 28 translates decompose K_{15x4} into T(6); the 183 blocks form an S(2,5,61)
    and agree with the block file shipped with the starter;
  * nauty: |Aut| of the design (point/block incidence graph) and pairwise non-isomorphism;
  * GF(2) rank of the 56 x 28 long-point crossing incidence matrix;
  * random admissible decorations (random Fano structure with quiet flags on a line,
    random phases with crossing product -1, random short labels, random special bijection):
      - degree-231 member is srg(496,231,102,112),
      - cell-constant kernel words (constant and character star) give srg(496,240,120,112),
      - a descendant is srg(495,238,109,119),
      - after a random Seidel switching and relabelling, the zero-sum 4-sets of rows of
        M = S - I are exactly the 980 affine planes inside cells.
Finally the Fano count (6 structures with a prescribed line) and the bound are evaluated.
"""
import sys
from fractions import Fraction
from itertools import combinations, permutations
from math import factorial, log10
from pathlib import Path
import random
import numpy as np
import pynauty

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'audit'))
sys.path.insert(0, str(HERE.parent / 'family'))
from independent_rebuild import (chi, part, tiles_from, read_starter, frame, crossing,
                                 regular_graph, described_z2, check_srg)
from fibre_planes import zero_sum_quads

STARTERS = HERE / 'starters'
LINES = [{a, b, a ^ b} for a, b in combinations(range(1, 8), 2) if a < (a ^ b) and b < (a ^ b)]
assert len(LINES) == 7


def blocks_of(tiles):
    B = [frozenset([4 * g, 4 * g + 1, 4 * g + 2, 4 * g + 3, 60]) for g in range(15)]
    for T in tiles:
        for a in range(6):
            B.append(frozenset(p for p, e in T.items() if a in e))
    assert len(B) == 183 and all(len(b) == 5 for b in B)
    cover = {}
    for i, b in enumerate(B):
        for pq in combinations(sorted(b), 2):
            assert pq not in cover, 'pair covered twice'
            cover[pq] = i
    assert len(cover) == 61 * 60 // 2, 'not an S(2,5,61)'
    return B


def read_blocks(path):
    return {frozenset(map(int, l.split())) for l in path.read_text().split('\n') if l.strip()}


def design_nauty(B):
    n = 61 + len(B)
    adj = {61 + i: sorted(b) for i, b in enumerate(B)}
    g = pynauty.Graph(n, adjacency_dict=adj, vertex_coloring=[set(range(61)), set(range(61, n))])
    _, size, exp, orbits, norb = pynauty.autgrp(g)
    order = round(size * 10 ** exp)
    inf_orbit = sum(1 for p in range(61) if orbits[p] == orbits[60])
    return pynauty.certificate(g), order, inf_orbit


def gf2_rank(rows):
    rows = [int(''.join(map(str, r)), 2) for r in rows]
    rank = 0
    for bit in reversed(range(max(r.bit_length() for r in rows))):
        piv = next((r for r in rows if (r >> bit) & 1), None)
        if piv is None:
            continue
        rows = [r ^ piv if (r >> bit) & 1 else r for r in rows if r != piv]
        rank += 1
    return rank


def crossing_matrix(tiles):
    M = []
    for p in range(56):
        M.append([1 if p in T and crossing(T[p]) else 0 for T in tiles])
        assert sum(M[-1]) == 4
    return M


def random_decoration(tiles, rng):
    chars, signs = [dict() for _ in range(60)], [dict() for _ in range(60)]
    for p in range(60):
        mine = [t for t in range(28) if p in tiles[t]]
        if p < 56:
            line = sorted(rng.choice(LINES))
            rest = [c for c in range(1, 8) if c not in line]
            rng.shuffle(line), rng.shuffle(rest)
            cross = [t for t in mine if crossing(tiles[t][p])]
            quiet = [t for t in mine if t not in cross]
            for t, c in zip(quiet + cross, line + rest):
                chars[p][t] = c
            s = [rng.choice((1, -1)) for _ in range(7)]
            for t, e in zip(mine, s):
                signs[p][t] = e
            prod = np.prod([signs[p][t] for t in cross])
            if prod == 1:
                signs[p][cross[0]] *= -1
        else:
            cs = list(range(1, 8))
            rng.shuffle(cs)
            for t, c in zip(mine, cs):
                chars[p][t], signs[p][t] = c, rng.choice((1, -1))
    return chars, signs


def special_permuted(F, rng):
    """Random bijection parts -> nonzero characters of F_2^4 on the special columns."""
    perm = list(range(1, 16))
    rng.shuffle(perm)
    F = F.copy()
    for j in range(15):
        F[j, :16] = [chi(perm[j], x) for x in range(16)]
    return F


def kernel_words(F, rng):
    """Constant-star and character-star words, constant on ordinary cells."""
    out = []
    for kind in ('constant', 'character'):
        w = np.ones(496, dtype=np.int64)
        active = rng.randrange(15) if kind == 'character' else None
        if active is not None:
            w[:16] = F[active, :16]
        for g in range(15):
            cells = [4 * g + i for i in range(4)]
            rng.shuffle(cells)
            neg = cells[:3] if g == active else cells[:2]
            for q in neg:
                w[16 + 8 * q:24 + 8 * q] = -1
        assert not np.any(F @ w), kind
        out.append((kind, w))
    return out


def check_fibres(S, rng):
    n = 496
    d = np.array([rng.choice((1, -1)) for _ in range(n)])
    perm = list(range(n))
    rng.shuffle(perm)
    S2 = (d[:, None] * S * d[None, :])[np.ix_(perm, perm)]
    A = ((S2 - np.eye(n, dtype=np.int64)) < 0).astype(np.int64) - np.eye(n, dtype=np.int64)  # M=-1 off diag
    Q = zero_sum_quads(A)
    cell = lambda v: -1 if v < 16 else (v - 16) // 8
    inside = sum(1 for q in Q if len({cell(perm[v]) for v in q}) == 1)
    return len(Q), inside


def main():
    rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 20261125)
    certs, ranks = [], []
    for name in ('design00', 'design01', 'design02'):
        tiles = tiles_from(read_starter(STARTERS / f'{name}.starter'))
        B = blocks_of(tiles)
        assert set(B) == read_blocks(STARTERS / f'{name}.blocks'), 'block file differs'
        cert, order, inf_orbit = design_nauty(B)
        certs.append(cert)
        r = gf2_rank(crossing_matrix(tiles))
        ranks.append(r)
        print(f'{name}: S(2,5,61) ok, block file agrees, |Aut(design)| = {order}, '
              f'orbit of the added point has size {inf_orbit}, crossing rank {r}')
        for trial in range(2):
            chars, signs = random_decoration(tiles, rng)
            F = special_permuted(frame(tiles, chars, signs), rng)
            A, w, S = regular_graph(F, described_z2())
            spec = check_srg(A)
            msg = [f'  member {trial}: 231 {spec}']
            for kind, kw in kernel_words(F, rng):
                assert np.array_equal(S @ kw, -15 * kw)
                B240 = ((kw[:, None] * S * kw[None, :]) > 0).astype(np.int64)
                check_srg(B240, 240, 120, 112)
                msg.append(f'240[{kind}] ok')
            r = rng.randrange(496)
            T = S * S[r][:, None] * S[r][None, :]
            keep = [v for v in range(496) if v != r]
            D = (T[np.ix_(keep, keep)] < 0).astype(np.int64)
            check_srg(D, 238, 109, 119)
            msg.append('495 descendant ok')
            if trial == 0:
                nq, inside = check_fibres(S, rng)
                msg.append(f'switched+relabelled zero-sum 4-sets {nq}, inside cells {inside}')
                assert nq == inside == 980
            print(', '.join(msg))
    assert len(set(certs)) == 3
    print('three designs pairwise non-isomorphic: True')

    # Fano structures on 7 labelled points: count, and those having a given 3-set as a line
    structures = set()
    for p in permutations(range(1, 8)):
        structures.add(frozenset(frozenset(p[c - 1] for c in L) for L in LINES))
    with_line = sum(1 for s in structures if frozenset({1, 2, 3}) in s)
    print(f'Fano structures on 7 points: {len(structures)}, with a prescribed line: {with_line}')
    assert len(structures) == 30 and with_line == 6
    assert factorial(15) % 20160 == 0

    base = Fraction(6 ** 56 * 30 ** 4 * (factorial(15) // 20160), 84)
    total = sum(base * 2 ** (364 - 184 - (28 - r)) for r in ranks)
    claimed = Fraction(6 ** 56 * 30 ** 4, 84) * Fraction(factorial(15), 20160) * 2 ** 152 * (2 ** 22 + 2 * 2 ** 19)
    assert total == claimed
    crude = sum(base * 2 ** (364 - 184 - 28) for _ in ranks)
    print(f'lower bound = {float(total):.4e} (log10 {log10(total.numerator) - log10(total.denominator):.3f}); '
          f'without the tile refinement {float(crude):.4e}')


if __name__ == '__main__':
    main()

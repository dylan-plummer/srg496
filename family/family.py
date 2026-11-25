"""Random members of the srg(496,231,102,112) family.

Free data, for the fixed starter and the cut {0,1}|{2,3,4,5}:
  long-sheet point p (56 of them), crossing in 4 tiles, quiet in 3:
    * a hyperplane H < F_2^3; quiet flags get H \\ {0}, crossing flags get F_2^3 \\ H
      (an affine plane avoiding 0), each in any order;
    * root signs: quiet flags arbitrary, crossing flags with product -1;
  short-sheet point (4 of them): any bijection tiles -> characters 1..7, any signs;
  special point: any bijection parts -> nonzero characters of F_2^4.
Every choice gives a tight frame (Walsh orthogonality per tile) and w = F^T z is a sign
vector, so every choice is an srg(496,231,102,112).  Usage:
    python family.py N [seed]      # N random members, nauty certificates and |Aut|
"""
import os
os.environ.setdefault('OMP_NUM_THREADS', '4')
import sys, random, hashlib
from pathlib import Path
import numpy as np
import pynauty

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'audit'))
from independent_rebuild import (tiles_from, read_starter, frame, crossing, chi, part,
                                 regular_graph, described_z2)

HYPERPLANES = [[x for x in range(1, 8) if bin(x & h).count('1') % 2 == 0] for h in range(1, 8)]

def random_labels(tiles, rng):
    chars, signs = [dict() for _ in range(60)], [dict() for _ in range(60)]
    for p in range(60):
        mine = [t for t in range(28) if p in tiles[t]]
        if p >= 56:
            cs = list(range(1, 8)); rng.shuffle(cs)
            for t, c in zip(mine, cs):
                chars[p][t], signs[p][t] = c, rng.choice((1, -1))
            continue
        cross = [t for t in mine if crossing(tiles[t][p])]
        quiet = [t for t in mine if not crossing(tiles[t][p])]
        H = rng.choice(HYPERPLANES)
        plane = [x for x in range(1, 8) if x not in H]
        rng.shuffle(plane); qs = list(H); rng.shuffle(qs)
        sg = [rng.choice((1, -1)) for _ in range(3)]
        sg.append(-sg[0] * sg[1] * sg[2])
        for t, c, s in zip(cross, plane, sg):
            chars[p][t], signs[p][t] = c, s
        for t, c in zip(quiet, qs):
            chars[p][t], signs[p][t] = c, rng.choice((1, -1))
    return chars, signs

def special_perm(F, rng):
    """Relabel the special columns: part j gets character pi(j)."""
    pi = list(range(1, 16)); rng.shuffle(pi)
    for x in range(16):
        for j in range(15):
            F[j, x] = chi(pi[j], x)
    return F

def member(tiles, rng):
    chars, signs = random_labels(tiles, rng)
    F = special_perm(frame(tiles, chars, signs), rng)
    A, w, S = regular_graph(F, described_z2())
    n = 496
    J = np.ones((n, n), dtype=np.int64)
    assert np.all(A.sum(1) == 231)
    assert np.array_equal(A @ A, 119 * np.eye(n, dtype=np.int64) - 10 * A + 112 * J)
    return A

def nauty(A):
    g = pynauty.Graph(len(A), adjacency_dict={i: list(np.nonzero(A[i])[0]) for i in range(len(A))})
    cert = hashlib.sha256(pynauty.certificate(g)).hexdigest()[:16]
    _, size, exp, _, norb = pynauty.autgrp(g)
    return cert, int(round(size * 10 ** exp)), norb

def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    rng = random.Random(seed)
    tiles = tiles_from(read_starter())
    certs = {}
    for i in range(N):
        A = member(tiles, rng)
        cert, aut, norb = nauty(A)
        certs.setdefault(cert, []).append(i)
        print(f'member {i}: srg(496,231,102,112) OK, cert {cert}, |Aut| {aut}, orbits {norb}', flush=True)
    print(f'{N} random members, {len(certs)} isomorphism classes')

if __name__ == '__main__':
    main()

"""Identify the regular two-graph of the 496 frame.

Compares |Aut| of the Taylor double cover (2|Aut(two-graph)|) with the symplectic
two-graph on the 496 singular vectors of an elliptic quadric in F_2^10, and
reports the descendant automorphism groups.
"""
import os
os.environ.setdefault('OMP_NUM_THREADS', '8')
import sys
from itertools import product
from pathlib import Path
import numpy as np
import pynauty

sys.path.insert(0, str(Path(__file__).resolve().parent))
from independent_rebuild import (tiles_from, read_starter, described_labels, frame)

def cover_graph(S):
    n = len(S)
    adj = {}
    for i in range(n):
        for e in (0, 1):
            nb = []
            for j in range(n):
                if j == i:
                    continue
                for d in (0, 1):
                    if (-1) ** (e + d) * S[i, j] == -1:
                        nb.append(2 * j + d)
            adj[2 * i + e] = nb
    return pynauty.Graph(2 * n, adjacency_dict=adj)

def group_order(g):
    _, size, exp, orbits, norb = pynauty.autgrp(g)
    return size * 10 ** exp, norb, orbits

def symplectic_seidel():
    # elliptic form on F_2^10: x1x2+...+x7x8 + x9^2+x9x10+x10^2
    def Q(v):
        s = sum(v[2*i] * v[2*i+1] for i in range(4)) + v[8] + v[8]*v[9] + v[9]
        return s % 2
    pts = [v for v in product((0, 1), repeat=10) if Q(v) == 0]
    assert len(pts) == 496
    P = np.array(pts)
    Bm = np.zeros((10, 10), dtype=int)
    for i in range(5):
        Bm[2*i, 2*i+1] = Bm[2*i+1, 2*i] = 1
    Bv = (P @ Bm @ P.T) % 2
    S = 1 - 2 * Bv
    np.fill_diagonal(S, 0)
    return S

def seidel_spectrum(S):
    ev = np.round(np.linalg.eigvalsh(S.astype(float))).astype(int)
    v, m = np.unique(ev, return_counts=True)
    return dict(zip(v.tolist(), m.tolist()))

def descendant(S, v):
    """Switch so that v is isolated, delete v: graph on the other 495 points."""
    s = S[v].copy(); s[v] = 1
    T = s[:, None] * S * s[None, :]
    keep = [i for i in range(len(S)) if i != v]
    T = T[np.ix_(keep, keep)]
    A = (T == -1)
    return pynauty.Graph(len(keep), adjacency_dict={i: list(np.nonzero(A[i])[0]) for i in range(len(keep))})

def main():
    """Descendant test.  The symplectic two-graph on 496 points is 2-transitive (Sp(10,2)), so all
    its descendants are isomorphic with |Aut| = |Sp(10,2)|/496.  Nauty on the 992-vertex double
    covers is slow, so we compare descendants instead (cover mode: python twograph.py cover)."""
    tiles = tiles_from(read_starter())
    F = frame(tiles, *described_labels(tiles))
    S = F.T @ F - 15 * np.eye(496, dtype=np.int64)
    Ssp = symplectic_seidel()
    print('frame two-graph Seidel spectrum', seidel_spectrum(S))
    print('symplectic two-graph Seidel spectrum', seidel_spectrum(Ssp), flush=True)
    if len(sys.argv) > 1 and sys.argv[1] == 'cover':
        for name, M in (('frame', S), ('symplectic', Ssp)):
            order, norb, orbits = group_order(cover_graph(M))
            print(f'{name}: |Aut(double cover)| = {order}, cover orbits {norb}', flush=True)
        return
    order, norb, _ = group_order(descendant(Ssp, 0))
    print(f'symplectic descendant at point 0: |Aut| = {order:.6g}, orbits {norb}', flush=True)
    for v in (0, 1, 15, 16, 17, 16 + 8 * 28, 16 + 8 * 56, 495):
        order, norb, _ = group_order(descendant(S, v))
        print(f'frame descendant at column {v}: |Aut| = {order}, orbits {norb}', flush=True)

if __name__ == '__main__':
    main()

"""Rebuild srg(496,231,102,112) from the written description only (PROOF.txt / note Sec. 2-3).

Checks, in exact integer arithmetic:
  * the 28 translates of the starter decompose K_{15x4} into copies of T(6);
  * F F^T = 48 I on the row space (F^T F has eigenvalues 48 and 0), rank 155;
  * S = F^T F - 15 I is a Seidel matrix;
  * w = F^T z is a +-1 vector with S w = 33 w;
  * A = (J - I - D S D)/2 satisfies A^2 = 119 I - 10 A + 112 J;
  * A equals the exported graph6 file.
"""
import sys
from itertools import combinations
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / 'data'
GRAPHS = HERE.parent / 'graphs'

def chi(j, x):
    return -1 if bin(j & x).count('1') % 2 else 1

def point(sheet, pos):
    """Point index: sheet 0 -> 0..27, sheet 1 -> 28..55, sheet 2 -> 56..59."""
    return 56 + pos if sheet == 2 else 28 * sheet + pos

def part(p):
    if p >= 56:
        return 14
    s, pos = divmod(p, 28)
    return 7 * s + pos // 4

def act(sheet, pos, t):
    """Translate by t = 4c + d in Z/7 x F_2^2."""
    c, d = divmod(t, 4)
    if sheet == 2:
        return pos ^ d
    a, b = divmod(pos, 4)
    return 4 * ((a + c) % 7) + (b ^ d)

def read_starter(path=DATA / 'starter.tsv'):
    rows = []
    for line in path.read_text().splitlines():
        if line.strip() and not line.startswith('#'):
            u, v, s, q = map(int, line.split())
            rows.append(((min(u, v), max(u, v)), s, q))
    assert sorted(e for e, _, _ in rows) == list(combinations(range(6), 2))
    return rows

def tiles_from(starter):
    tiles = []
    for t in range(28):
        tiles.append({point(s, act(s, q, t)): e for e, s, q in starter})
    # decomposition of K_{15x4}
    seen = {}
    for t, T in enumerate(tiles):
        assert len(T) == 15
        for p, q in combinations(T, 2):
            if set(T[p]) & set(T[q]):
                key = (min(p, q), max(p, q))
                assert key not in seen
                seen[key] = t
    want = {(p, q) for p, q in combinations(range(60), 2) if part(p) != part(q)}
    assert set(seen) == want, 'not a decomposition'
    return tiles

def frame(tiles, chars=None, signs=None):
    """chars[p][t], signs[p][t]: character in 1..7 and sign of the root for point p in tile t."""
    F = np.zeros((15 + 6 * 28, 496), dtype=np.int64)
    for x in range(16):
        for j in range(15):
            F[j, x] = chi(j + 1, x)
    for p in range(60):
        base = 16 + 8 * p
        mine = [t for t in range(28) if p in tiles[t]]
        assert len(mine) == 7
        assert sorted(chars[p][t] for t in mine) == list(range(1, 8))
        for x in range(8):
            col = base + x
            F[part(p), col] = 1
            for t in mine:
                a, b = tiles[t][p]
                v = signs[p][t] * chi(chars[p][t], x)
                F[15 + 6 * t + a, col] += v
                F[15 + 6 * t + b, col] -= v
    return F

def crossing(e):
    a, b = e
    return (a < 2) != (b < 2)

def described_labels(tiles):
    """Section 3: crossing flags get 4,5,6,7, quiet flags 1,2,3, each in tile order;
    root e_a-e_b except the crossing flag with character 7; short sheet 1..7 in tile order."""
    chars, signs = [dict() for _ in range(60)], [dict() for _ in range(60)]
    for p in range(60):
        mine = [t for t in range(28) if p in tiles[t]]
        cross = [t for t in mine if crossing(tiles[t][p])]
        quiet = [t for t in mine if not crossing(tiles[t][p])]
        if p < 56:
            assert len(cross) == 4
            for c, t in zip((4, 5, 6, 7), cross):
                chars[p][t], signs[p][t] = c, (-1 if c == 7 else 1)
            for c, t in zip((1, 2, 3), quiet):
                chars[p][t], signs[p][t] = c, 1
        else:
            assert not cross
            for c, t in zip(range(1, 8), mine):
                chars[p][t], signs[p][t] = c, 1
    return chars, signs

def regular_graph(F, z2):
    """z2 = 2z. Returns (A, w, S) or raises."""
    gram = F.T @ F
    assert np.array_equal(gram @ gram, 48 * gram)
    assert np.trace(gram) == 48 * 155
    S = gram - 15 * np.eye(496, dtype=np.int64)
    off = S[~np.eye(496, dtype=bool)]
    assert np.all(np.abs(off) == 1) and np.all(np.diag(S) == 0)
    w2 = F.T @ z2
    assert np.all(np.abs(w2) == 2), 'w is not a sign vector'
    w = w2 // 2
    assert np.array_equal(S @ w, 33 * w)
    T = w[:, None] * S * w[None, :]
    A = ((1 - np.eye(496, dtype=np.int64)) - T) // 2
    return A, w, S

def described_z2():
    z2 = np.zeros(15 + 6 * 28, dtype=np.int64)
    z2[14] = 2
    for t in range(28):
        z2[15 + 6 * t] = z2[15 + 6 * t + 1] = 1
    return z2

def check_srg(A, k=231, lam=102, mu=112):
    n = len(A)
    assert np.array_equal(A, A.T) and np.all(np.diag(A) == 0) and set(np.unique(A)) <= {0, 1}
    assert np.all(A.sum(1) == k)
    J = np.ones((n, n), dtype=np.int64)
    assert np.array_equal(A @ A, (k - mu) * np.eye(n, dtype=np.int64) + (lam - mu) * A + mu * J)
    ev = np.round(np.linalg.eigvalsh(A.astype(float))).astype(int)
    vals, mult = np.unique(ev, return_counts=True)
    return dict(zip(vals.tolist(), mult.tolist()))

def read_g6(path):
    data = path.read_bytes().strip()
    n = data[0] - 63
    assert n == 63
    n = ((data[1] - 63) << 12) | ((data[2] - 63) << 6) | (data[3] - 63)
    bits = []
    for c in data[4:]:
        v = c - 63
        bits.extend((v >> s) & 1 for s in range(5, -1, -1))
    A = np.zeros((n, n), dtype=np.int64)
    k = 0
    for j in range(n):
        for i in range(j):
            A[i, j] = A[j, i] = bits[k]
            k += 1
    return A

def main():
    tiles = tiles_from(read_starter())
    chars, signs = described_labels(tiles)
    F = frame(tiles, chars, signs)
    A, w, S = regular_graph(F, described_z2())
    spec = check_srg(A)
    print('srg(496,231,102,112) from the description; spectrum', spec)
    ev = np.round(np.linalg.eigvalsh(S.astype(float))).astype(int)
    print('Seidel spectrum', dict(zip(*[x.tolist() for x in np.unique(ev, return_counts=True)])))
    B = read_g6(GRAPHS / 'srg_496_231_102_112.g6')
    print('equal to exported graph6:', np.array_equal(A, B))
    C = read_g6(GRAPHS / 'srg_496_264_144_136.g6')
    print('complement file is J-I-A:', np.array_equal(C, 1 - np.eye(496, dtype=np.int64) - A))
    print('complement spectrum', check_srg(C, 264, 144, 136))

if __name__ == '__main__':
    main()

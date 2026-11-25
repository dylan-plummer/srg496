"""Translation-equivariant Walsh decorations of the 496-column frame."""
from itertools import permutations
from pathlib import Path
import sys

import networkx as nx
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'seed'))
import construct as base


def mul8(x, y):
    """Multiply modulo X^3 + X + 1, representing coefficients by bits."""
    result = 0
    while y:
        if y & 1:
            result ^= x
        y //= 2
        x <<= 1
        if x & 8:
            x ^= 11
    return result


def power(i):
    value = 1
    for _ in range(i % 7):
        value = mul8(value, 2)
    return value


def trace(x):
    square = mul8(x, x)
    return x ^ square ^ mul8(square, square)


def chi(c, x):
    return 1 - 2 * trace(mul8(c, x))


def rotate(x):
    return 4 * (2 * (x // 4) % 7) + (0, 2, 3, 1)[x % 4]


def minus(x, y):
    return 4 * ((x // 4 - y // 4) % 7) + ((x ^ y) & 3)


def palettes(tiles, sheet):
    origin = 28 * sheet
    incident = [t for t, tile in enumerate(tiles) if origin in tile]
    crossing = {t for t in incident if min(tiles[t][origin]) < 2 <= max(tiles[t][origin])}
    options = []
    for assignment in permutations(range(1, 8)):
        labels = dict(zip(incident, assignment))
        if all(trace(c) == (t in crossing) and labels[rotate(t)] == mul8(c, c)
               for t, c in labels.items()):
            options.append(labels)
    return options


def make(first=0, second=0):
    tiles = base.cover()
    labels = [palettes(tiles, s)[i] for s, i in enumerate((first, second))]
    F = np.zeros((183, 496), dtype=np.int64)
    for sheet in range(2):
        for i in range(7):
            for x in range(16):
                F[7 * sheet + i, x] = chi(power(i), x % 8) * (1 - 2 * (x // 8) * sheet)
    F[14, :16] = [1 - 2 * (x // 8) for x in range(16)]

    for p in range(60):
        sheet = min(p // 28, 2)
        origin = p % 28 if sheet < 2 else p - 56
        cols = slice(16 + 8 * p, 24 + 8 * p)
        F[p // 4, cols] = 1
        for t, tile in enumerate(tiles):
            if p not in tile:
                continue
            a, b = tile[p]
            # The moving K6 symbols carry the cyclic orientation 3,4,5.
            if a >= 3 and (b - a) % 3 != 1:
                a, b = b, a
            c = labels[sheet][minus(t, origin)] if sheet < 2 else power(t // 4)
            phase = -1 if sheet < 2 and c == 1 else 1
            values = [phase * chi(c, x) for x in range(8)]
            start = 15 + 6 * t
            F[start + a, cols] = values
            F[start + b, cols] = -np.asarray(values)
    return F


def switch(F, cut=(0, 1), special=14):
    z = np.zeros(183, dtype=np.int64)
    z[special] = 2
    for t in range(28):
        for a in cut:
            z[15 + 6 * t + a] = 1
    word = F.T @ z
    assert set(word) == {-2, 2}
    w = word // 2
    S = F.T @ F - 15 * np.eye(496, dtype=np.int64)
    A = (w[:, None] * S * w[None, :] < 0).astype(np.int64)
    base.verify(A, 231, 102, 112)
    return A, S, w


if __name__ == '__main__':
    A, _, _ = switch(make())
    nx.write_graph6(nx.from_numpy_array(A), HERE / 'srg496_symmetric.g6', header=False)

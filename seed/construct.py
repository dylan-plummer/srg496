"""Construct the two 496-vertex graphs and their supporting 2-(61,5,1) design."""
from collections import Counter
from itertools import combinations
from pathlib import Path

import networkx as nx
import numpy as np

HERE = Path(__file__).resolve().parent
GROUPS = 15
TILES = 28
SPECIAL = 60

def group(p):
    return p // 4

def translate(sheet, pos, t):
    bit_part = (pos ^ t) & 3
    if sheet == 2:
        return 56 + bit_part
    return sheet * 28 + 4 * ((pos // 4 + t // 4) % 7) + bit_part

def walsh(a, b):
    return 1 if (a & b).bit_count() % 2 == 0 else -1

def cover():
    starter = [tuple(map(int, s.split())) for s in (HERE / 'starter.tsv').read_text().splitlines()
               if s and not s.startswith('#')]
    assert {tuple(sorted((a,b))) for a,b,_,_ in starter} == set(combinations(range(6),2))
    tiles = []
    count = Counter()
    for t in range(TILES):
        tile = {translate(sheet, pos, t): (a, b) for a, b, sheet, pos in starter}
        assert len(tile) == 15
        for p,q in combinations(sorted(tile),2):
            if set(tile[p]) & set(tile[q]):
                count[p,q] += 1
        tiles.append(tile)
    assert count == Counter({(p,q):1 for p,q in combinations(range(60),2) if group(p) != group(q)})
    return tiles

def design(tiles):
    blocks = [tuple([*range(4*j,4*j+4), SPECIAL]) for j in range(GROUPS)]
    for tile in tiles:
        for u in range(6):
            blocks.append(tuple(sorted(p for p,e in tile.items() if u in e)))
    assert len(set(blocks)) == 183 and all(len(b)==5 for b in blocks)
    assert Counter(e for b in blocks for e in combinations(b,2)) == Counter(combinations(range(61),2))
    return blocks

def frame(tiles):
    # Rows are the 15 groups, followed by six coordinates in each tile.
    F = np.zeros((GROUPS+6*TILES, 496), dtype=np.int64)
    for x in range(16):
        F[:GROUPS, x] = [walsh(j+1,x) for j in range(GROUPS)]
    for p in range(60):
        cols = slice(16+8*p, 24+8*p)
        F[group(p),cols] = 1
        flags = [(t,tile[p]) for t,tile in enumerate(tiles) if p in tile]
        active = [(t,e) for t,e in flags if (e[0]<2) != (e[1]<2)]
        quiet = [(t,e) for t,e in flags if (t,e) not in active]
        assert len(flags)==7 and len(active)==(4 if p<56 else 0)
        for j,(t,(a,b)) in enumerate(quiet+active,1):
            phase = 1
            if (a<2) != (b<2):
                phase = (1 if a<2 else -1) * (-1 if j==7 else 1)
            row = GROUPS + 6*t
            values = [phase*walsh(j,x) for x in range(8)]
            F[row+a,cols] = values
            F[row+b,cols] = [-c for c in values]
    return F

def verify(A, k, lam, mu):
    n = len(A)
    assert set(np.unique(A)) == {0,1} and np.all(A.diagonal()==0) and np.array_equal(A,A.T)
    assert np.all(A.sum(axis=1)==k)
    assert np.array_equal(A@A, (k-mu)*np.eye(n,dtype=np.int64) + (lam-mu)*A + mu)

def main():
    tiles = cover()
    blocks = design(tiles)
    F = frame(tiles)
    gram = F.T @ F
    assert np.array_equal(gram @ gram,48*gram)
    S = gram - 15*np.eye(496,dtype=np.int64)
    assert np.array_equal(S.diagonal(),np.zeros(496)) and np.all(abs(S + np.eye(496, dtype=int))==1)
    # Use twice the prescribed row-space vector, so all arithmetic is integral.
    z2 = np.zeros(len(F),dtype=np.int64)
    z2[14]=2
    for t in range(TILES):
        z2[GROUPS+6*t:GROUPS+6*t+2] = 1
    w2 = F.T @ z2
    assert set(w2)=={-2,2}
    w = w2//2
    assert np.array_equal(S @ w,33*w)
    T = w[:,None]*S*w[None,:]
    A = (1-np.eye(496,dtype=np.int64)-T)//2
    C = 1-np.eye(496,dtype=np.int64)-A
    for M,k,lam,mu in [(A,231,102,112),(C,264,144,136)]:
        verify(M,k,lam,mu)
        nx.write_graph6(nx.from_numpy_array(M), HERE / f'srg_496_{k}_{lam}_{mu}.g6',header=False)
        print(f'({len(M)},{k},{lam},{mu}): {M.sum()//2} edges, {len(M)*(len(M)-1)//2} pairs verified')
    with (HERE/'steiner_2_61_5.txt').open('w') as out:
        for b in blocks:
            print(*b, file=out)
    print(f'{len(blocks)} Steiner blocks, {len(tiles)} tiles, frame rank {gram.trace()//48}')
if __name__ == '__main__':
    main()

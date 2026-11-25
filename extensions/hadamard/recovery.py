"""Check the starter-only conditions for recognizing the 976-point constructions."""
from collections import Counter
from itertools import combinations,permutations
from math import factorial
from pathlib import Path
import json
from fractions import Fraction

from build import HERE, START, tiles


def full_recovery(design):
    blocks=[set(map(int,l.split())) for l in (START/f'design{design:02}.blocks').read_text().splitlines()]
    ts=tiles(design)
    split_witness=[]
    for p,q in combinations(range(61),2):
        left=[u-{p} for u in blocks if p in u and q not in u]
        right=[u-{q} for u in blocks if q in u and p not in u]
        neighbours=[{j for j,v in enumerate(right) if u&v} for u in left]
        assert len(left)==len(right)==14 and all(len(n)==4 for n in neighbours)
        reached={0}
        while True:
            through=set.union(*(neighbours[i] for i in reached))
            again={i for i,n in enumerate(neighbours) if n&through}
            if reached==again:
                break
            reached=again
        assert len(reached)==14
        if q<60 and p//4==q//4:
            def port(x,r):
                hits=[j for j,t in enumerate(ts) if x in t and r in t and set(t[x])&set(t[r])]
                only,=hits
                return only
            split_witness.append(sum(port(p,r)!=port(q,r) for r in range(60) if r//4!=p//4))
    assert len(split_witness)==90 and min(split_witness)>0
    print('design',design,'all 1830 comparison graphs connected; witnesses for 90 same-part pairs',dict(Counter(split_witness)))


def main():
    for d in range(3):
        full_recovery(d)
    structures={}
    for perm in permutations(range(1,8)):
        fano=tuple((i,j,k) for i,j,k in combinations(range(7),3) if perm[i]^perm[j]^perm[k]==0)
        words={tuple((c&x).bit_count()%2 for c in perm) for x in range(8)}
        structures[fano]=words
    base=next(iter(structures.values()))
    print('number of choices by size of intersection',dict(Counter(len(w&base) for w in structures.values())))
    assert len(structures)==30
    print('distinct parents of one fixed half, lower bound',(30**60+83)//84)
    rec=json.loads((HERE/'hadamards.json').read_text())
    lower=Fraction(6**56*30**4*factorial(15)*2**(364-180)*(2**22+2*2**19),84*2**28)
    old=lower/rec[0]['order']
    gains=[lower/r['order'] for r in rec]
    print('relative weights',[v/old for v in gains])
    print('old rational bound',old,'new four-type rational bound',sum(gains[1:]),'combined',sum(gains))
    print('scientific, new and combined',float(sum(gains[1:])),float(sum(gains)))
    assert sum(gains[1:])==485*old


if __name__=='__main__':
    main()

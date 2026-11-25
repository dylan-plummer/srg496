"""Enumerate regular switches arising from cuts of the six-symbol tiles."""
from collections import Counter
from functools import reduce
from itertools import combinations, product
from operator import xor
import csv
import gzip
import numpy as np
from family import HERE, base, make, switch, rotate
from cuts import setup, nullspace


def linear_solutions(equations, n):
    """Affine solutions to parity rows (mask, right-hand side), as bitsets."""
    pivots={}
    for v,r in equations:
        while v:
            k=v.bit_length()-1
            if k not in pivots:
                pivots[k]=v,r
                break
            prior,b=pivots[k]
            v ^= prior
            r ^= b
        if not v and r:
            return []
    seeds=[0]
    for j in range(n):
        if j not in pivots:
            seeds += [a | (1<<j) for a in seeds]
    for k,(v,r) in sorted(pivots.items()):
        seeds=[a | ((((v&a).bit_count() ^ r)&1)<<k) for a in seeds]
    return seeds


def all_cuts(flags, F, group):
    equations=[]
    for p,row in enumerate(flags):
        es=dict(row)
        for c,e in es.items():
            rhs=0 if p//4==group else reduce(xor,(es[i] for i in (1,2,4) if c&i),0)
            equations.append((e ^ rhs,0))
    tiles=base.cover()
    raw={}
    for t,T in enumerate(tiles):
        for p,(a,b) in T.items():
            # Entry at the zero Walsh column includes the prescribed root orientation.
            raw[t,p]=(a,b,int(F[15+6*t+a,16+8*p]))
    for cut in linear_solutions(equations,140):
        bits=[(cut>>(5*t))&31 for t in range(28)]
        sign_rows=[(1<<t,0) for t,b in enumerate(bits) if not b]
        for p in range(60):
            if p//4==group:
                continue
            mask=0
            parity=1
            for t,T in enumerate(tiles):
                if p not in T:
                    continue
                a,b,orientation=raw[t,p]
                difference=((bits[t]<<1)>>a&1)-((bits[t]<<1)>>b&1)
                if difference:
                    mask |= 1<<t
                    parity ^= difference*orientation < 0
            if mask.bit_count()!=4:
                break
            sign_rows.append((mask,parity))
        else:
            for orientations in linear_solutions(sign_rows,28):
                yield (group,bytes(b+32*((orientations>>t)&1) for t,b in enumerate(bits)))


def group_elements():
    rotbits=(0,2,3,1)
    for r,b,k in product(range(7),range(4),range(3)):
        def act(t):
            for _ in range(k):
                t=rotate(t)
            return 4*((t//4+r)%7)+(t%4 ^ b)
        tiles=[act(t) for t in range(28)]
        groups=[7*s+act(4*i)//4 for s in range(2) for i in range(7)]+[14]
        pair=[x^b for x in range(4)]
        for j in range(4):
            y=j
            for _ in range(k):
                y=rotbits[y]
            pair[j]=y^b
        lut=[]
        for value in range(64):
            new=value&32
            for s in range(1,6):
                q=3+(s-3+k)%3 if s>=3 else s
                new |= ((value>>(s-1))&1)<<(q-1)
            lut.append(new)
        yield (r,b,k),tiles,groups,pair,lut


def transform(key,action):
    g,bits=key
    _,tiles,groups,_,lut=action
    new=bytearray(28)
    for t,y in enumerate(bits):
        new[tiles[t]]=lut[y]
    return groups[g],bytes(new)


def balanced_count():
    ordinary=0
    pointed=0
    for _,_,groups,perm,_ in group_elements():
        seen=set()
        factors=[]
        for g in range(15):
            if g in seen:
                continue
            h=g
            p=list(range(4))
            while h not in seen:
                seen.add(h)
                p=[perm[q] for q in p]
                h=groups[h]
            pairs=sum({p[a],p[b]}=={a,b} for a,b in combinations(range(4),2))
            singleton=sum(p[a]==a for a in range(4)) if groups[g]==g else 0
            factors.append((pairs,singleton))
        from math import prod
        ordinary+=prod(a for a,b in factors)
        pointed+=sum(b*prod(c for j,(c,d) in enumerate(factors) if j!=i) for i,(a,b) in enumerate(factors))
    assert ordinary%84==pointed%84==0
    return ordinary//84,pointed//84


def graph_from_code(F,code):
    g,bits=code
    z2=np.zeros(183,dtype=np.int64)
    z2[g]=2
    for t,b in enumerate(bits):
        for a in range(1,6):
            z2[15+6*t+a]=(1-2*(b>>5))*((b>>(a-1))&1)
    w2=F.T@z2
    assert set(w2)=={-2,2}
    S=F.T@F-15*np.eye(496,dtype=np.int64)
    return ((w2[:,None]*S*w2[None,:])<0).astype(np.int64)


def main():
    actions=list(group_elements())
    counts=Counter()
    print('balanced kernel switches per two-graph',balanced_count(),flush=True)
    with gzip.open(HERE/'cut_classes.tsv.gz','wt') as out:
        print('palette\tgroup\tcut_hex\tautomorphisms',file=out)
        for i,j in product(range(3),repeat=2):
            F=make(i,j)
            flags=setup(i,j)
            whole=set()
            for g in (0,7,14):
                for x in all_cuts(flags,F,g):
                    whole.update(transform(x,a) for a in actions)
            hist=Counter()
            total=len(whole)
            while whole:
                x=min(whole)
                orbit={transform(x,a) for a in actions}
                assert orbit<=whole
                whole.difference_update(orbit)
                aut=84//len(orbit)
                hist[aut]+=1
                print(f'{i}{j}',x[0],x[1].hex(),aut,sep='\t',file=out)
            print(i,j,total,sum(hist.values()),sorted(hist.items()),flush=True)
            counts.update(hist)
    print('total',sum(counts.values()),sorted(counts.items()))


if __name__=='__main__':
    main()

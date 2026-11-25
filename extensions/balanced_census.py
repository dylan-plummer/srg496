"""Verify the signed-frame symmetry and both types of fibre-constant kernel switch."""
from collections import Counter
from functools import reduce
from itertools import combinations, product
from math import prod
import networkx as nx
import numpy as np
from family import HERE,base,make,switch,power,mul8
from switch_classes import group_elements,balanced_count


def permutation(element):
    (r,b,k),tiles,groups,colour,_=element
    p=list(range(496))
    for v in range(496):
        if v<16:
            x=v%8
            for _ in range(k):
                x=mul8(x,x)
            p[v]=v//8*8+mul8(power(-r),x)
        else:
            q,x=divmod(v-16,8)
            for _ in range(k):
                x=mul8(x,x)
            if q<56:
                q=(q//28)*28+tiles[q%28]
            else:
                q=56+colour[q-56]
                x=mul8(power(-r),x)
            p[v]=16+8*q+x
    return tuple(p)


def count_fixed(P):
    """Read the action on individual column fibres, directly from a vertex permutation."""
    cells=[(P[16+8*q]-16)//8 for q in range(60)]
    cycles=[]
    seen=set()
    for g in range(15):
        if g in seen:
            continue
        point=list(range(4*g,4*g+4))
        groups=0
        while point[0]//4 not in seen:
            seen.add(point[0]//4)
            point=[cells[x] for x in point]
            groups+=1
        two=sum({point[a],point[b]}=={4*g+a,4*g+b} for a,b in combinations(range(4),2))
        one=sum(point[a]==4*g+a for a in range(4)) if groups==1 else 0
        cycles.append((two,one))
    star=prod(x for x,y in cycles)
    character=sum(y*prod(u for j,(u,v) in enumerate(cycles) if j!=i) for i,(x,y) in enumerate(cycles))
    return star,character


def main():
    group=list(group_elements())
    generators=[permutation(a) for a in group if a[0] in ((1,0,0),(0,1,0),(0,2,0),(0,0,1))]
    perms={permutation(a) for a in group}
    obtained={tuple(range(496))}
    queue=list(obtained)
    for p in queue:
        for g in generators:
            q=tuple(g[i] for i in p)
            if q not in obtained:
                queue.append(q)
                obtained.add(q)
    assert len(perms)==84 and obtained==perms
    fixed=Counter(count_fixed(p) for p in perms)
    N=tuple(sum(k[i]*f for k,f in fixed.items())//84 for i in (0,1))
    assert N==balanced_count()
    print('Burnside fixed-word counts (constant star, character star), frequencies:',sorted(fixed.items()))
    print('per two-graph',N,'sum',sum(N),'all nine',9*sum(N))
    for i,j in product(range(3),repeat=2):
        F=make(i,j)
        A,S,_=switch(F)
        for p in generators:
            assert np.array_equal(S[np.ix_(p,p)],S) and np.array_equal(A[np.ix_(p,p)],A)
        if i==j==0:
            for active in range(15):
                w=np.ones(496,dtype=np.int64)
                w[:16]=F[active,:16]
                for g in range(15):
                    start=16+32*g
                    w[start:start+(24 if g==active else 16)]=-1
                assert np.all(F@w==0)
                if active==0:
                    B=(S*w[:,None]*w[None,:]>0).astype(np.int64)
                    base.verify(B,240,120,112)
                    nx.write_graph6(nx.from_numpy_array(B),HERE/'graphs'/'496_240_character.g6',header=False)
    print('All four generators preserve all nine raw matrices and degree-231 graphs.')
    print('All 15 character choices have been checked on the first seed.')


if __name__=='__main__':
    main()

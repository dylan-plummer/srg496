"""Solve the linear constraints on binary symbol cuts for a one-character star."""
from itertools import product
from family import base, palettes, minus, power

def nullspace(equations,n):
    pivots={}
    for v in equations:
        while v:
            k=v.bit_length()-1
            if k in pivots:
                v ^= pivots[k]
            else:
                pivots[k]=v
                break
    basis=[]
    for j in range(n):
        if j in pivots:
            continue
        v=1<<j
        for k,eq in sorted(pivots.items()):
            if (eq&v).bit_count()%2:
                v |= 1<<k
        basis.append(v)
    return basis

def setup(i,j):
    tiles=base.cover()
    pal=[palettes(tiles,s)[a] for s,a in enumerate((i,j))]
    columns=[[] for _ in range(60)]
    for t,T in enumerate(tiles):
        for p,(a,b) in T.items():
            c=pal[p//28][minus(t,p%28)] if p<56 else power(t//4)
            edge=sum(1<<(5*t+k-1) for k in (a,b) if k)
            columns[p].append((c,edge))
    return columns

def scan(columns,group):
    eq=[]
    for p,flags in enumerate(columns):
        es={c:e for c,e in flags}
        for c,e in es.items():
            if group==p//4:
                eq.append(e)
            else:
                rhs=0
                for i in (1,2,4):
                    if c&i: rhs ^= es[i]
                eq.append(e ^ rhs)
    basis=nullspace(eq,140)
    count=0
    if len(basis)<20:
        vs=[0]
        for b in basis:
            vs += [v^b for v in vs]
        for v in vs:
            if all(any((v&e).bit_count()%2 for c,e in flags) for p,flags in enumerate(columns) if p//4!=group):
                count+=1
    return len(basis),count

if __name__=='__main__':
    for i,j in product(range(3),repeat=2):
        cols=setup(i,j)
        print(i,j,[scan(cols,g) for g in (0,7,14)],flush=True)

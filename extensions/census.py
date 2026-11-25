"""Identify distinct regular switches using canonical labels on uncoloured graphs."""
from collections import Counter
from hashlib import sha256
from itertools import product
import time
from pathlib import Path
import sys
import numpy as np
import networkx as nx
import igraph as ig
from family import HERE,make,switch
def canonical(A):
    G=ig.Graph(n=len(A),edges=list(zip(*np.nonzero(np.triu(A,1)))))
    p=G.canonical_permutation()
    v=np.asarray(p)
    edges=np.packbits(A[np.ix_(v,v)][np.triu_indices(len(A),1)])
    return sha256(edges.tobytes()).hexdigest(),int(G.count_automorphisms())
if __name__=='__main__':
    seen={}; start=time.monotonic()
    with (HERE/'palettes.tsv').open('w') as f:
        print('sheet0\tsheet1\tcanonical_sha256\tgroup_order',file=f)
        for i,j in product(range(9),repeat=2):
            F=make(i,j);A,_,_=switch(F);h,g=canonical(A)
            print(i,j,h,g,sep='\t',file=f,flush=True)
            seen.setdefault(h,g)
    print('Distinct',len(seen),'full group orders',dict(Counter(seen.values())),'seconds',time.monotonic()-start)

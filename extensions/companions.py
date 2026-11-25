"""Classify descendants and count balanced kernel switches in the nine examples."""
from collections import Counter, defaultdict
import csv
from hashlib import sha256
from itertools import combinations
from math import ceil
import networkx as nx
import numpy as np
import pynauty
from family import HERE, make, switch
from census import canonical


def nauty_graph(A):
    return pynauty.Graph(len(A), adjacency_dict={i: np.flatnonzero(r).tolist() for i,r in enumerate(A)})


def descendant(S, v):
    signs = S[v].copy()
    signs[v] = 1
    T = signs[:,None] * S * signs[None,:]
    idx = np.arange(len(S)) != v
    return (T[np.ix_(idx,idx)] < 0).astype(np.int64)


def balanced(F,S, choices):
    w = np.ones(496, dtype=np.int64)
    for g, pair in enumerate(choices):
        for i in pair:
            lo = 16 + (4*g+i)*8
            w[lo:lo+8] = -1
    assert np.all(F @ w == 0)
    T = w[:,None] * S * w[None,:]
    return (T > 0).astype(np.int64)


def main():
    grid = list(csv.DictReader((HERE/'palettes.tsv').read_text().splitlines(), delimiter='\t'))
    classes = {}
    roots = {}
    exported = {}
    (HERE/'graphs').mkdir(exist_ok=True)
    with (HERE/'descendants.tsv').open('w') as desc, (HERE/'classes.tsv').open('w') as summary:
        print('palette\troot\tmultiplicity\tbliss_sha256\tpoint_stabilizer',file=desc)
        print('palette\tnauty_sha256\tgraph_group\ttwo_graph_group\troot_orbits\ttwo_graph_key',file=summary)
        for row in grid:
            h = row['canonical_sha256']
            if h in classes:
                continue
            i,j = map(int,(row['sheet0'],row['sheet1']))
            name = f'{i}{j}'
            F=make(i,j)
            A,S,_=switch(F)
            G=nauty_graph(A)
            _,order,e,orbits,_=pynauty.autgrp(G)
            assert int(order*10**e)==84
            cert=sha256(pynauty.certificate(G)).hexdigest()
            assert cert not in classes.values()
            classes[h]=cert
            nx.write_graph6(nx.from_numpy_array(A),HERE/'graphs'/f'496_231_{name}.g6',header=False)
            by_root=defaultdict(list)
            for root, mult in sorted(Counter(orbits).items()):
                D = descendant(S,root)
                key,stab=canonical(D)
                by_root[key].append((root,mult,stab))
            orders=set()
            for key,values in by_root.items():
                root=values[0][0]
                stab=values[0][2]
                assert all(v[2]==stab for v in values)
                size=sum(v[1] for v in values)
                orders.add(size*stab)
                print(name,root,size,key,stab,sep='\t',file=desc,flush=True)
                roots.setdefault(key,name)
            assert len(orders)==1
            group=orders.pop()
            sig=min(by_root)
            print(name,cert,84,group,len(by_root),sig,sep='\t',file=summary,flush=True)
            if sig not in exported:
                exported[sig]=group
                B = balanced(F,S,[(0,1)]*15)
                nx.write_graph6(nx.from_numpy_array(B),HERE/'graphs'/f'496_240_{name}.g6',header=False)
                nx.write_graph6(nx.from_numpy_array(descendant(S,0)),HERE/'graphs'/f'495_238_{name}.g6',header=False)
            print(name,'two-graph group',group,'descendants',len(by_root),flush=True)
    print('two-graphs', len(exported), 'distinct descendants',len(roots),'balanced switches lower bound',sum(ceil(6**15/order) for order in exported.values()),flush=True)


if __name__=='__main__':
    main()

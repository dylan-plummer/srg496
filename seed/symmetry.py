"""Compute the automorphism group of the uncolored graph with nauty and bliss."""
from pathlib import Path
import igraph as ig
import networkx as nx
import pynauty
from collections import Counter

src=Path(__file__).resolve().parents[1]/'graphs'/'srg_496_231_102_112.g6'
g=nx.read_graph6(src)
gens, size, exponent, orbit, n_orbits=pynauty.autgrp(
    pynauty.Graph(len(g),adjacency_dict={i:list(g[i]) for i in g}))
ident=tuple(range(len(g)))
elts={ident}
pending=[ident]
for a in pending:
    for b in gens:
        image=tuple(b[i] for i in a)
        if image not in elts:
            elts.add(image)
            pending.append(image)
assert len(elts)==int(size*10**exponent)==ig.Graph(n=len(g),edges=list(g.edges())).count_automorphisms()
orders=[]
for p in elts:
    assert all(g.has_edge(p[i],p[j]) for i,j in g.edges())
    q=p;order=1
    while q!=ident:
        q=tuple(p[j] for j in q)
        order+=1
    orders.append(order)
print('Full group order',len(elts),'element orders',sorted(orders))
print('Vertex orbit sizes and frequencies',sorted(Counter(Counter(orbit).values()).items()))

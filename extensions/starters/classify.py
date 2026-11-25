"""Canonicalize the Steiner incidence graphs of the normalized difference starters."""
from collections import Counter
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from family import base
from companions import nauty_graph
from export import graph6
import networkx as nx
import numpy as np
import pynauty


def tiles_from(rows):
    tiles=[]
    for t in range(28):
        tile={base.translate(s,x,t):(a,b) for a,b,s,x in rows}
        assert len(tile)==15
        tiles.append(tile)
    return tiles


def design_graph(blocks):
    adj={i:[] for i in range(244)}
    for i,b in enumerate(blocks,61):
        for p in b:
            adj[i].append(p)
            adj[p].append(i)
    return pynauty.Graph(244,adjacency_dict=adj)


def canonical(g):
    return sha256(pynauty.certificate(g)).hexdigest()


def group_order(g):
    _,a,b,_,_=pynauty.autgrp(g)
    return int(a*10**b)


def main(data):
    counts=Counter()
    designs={}
    original=canonical(design_graph(base.design(base.cover())))
    with (HERE/'class_index.tsv').open('w') as out:
        print('line\taction\tdesign_nauty_sha256',file=out)
        for i,line in enumerate(data.read_text().splitlines(),1):
            mode,*flat=map(int,line.split())
            rows=[]
            for mask,s,x in zip(flat[::3],flat[1::3],flat[2::3]):
                a,b=(v for v in range(6) if mask>>v&1)
                rows.append((a,b,s,x))
            tiles=tiles_from(rows)
            blocks=base.design(tiles)
            design=design_graph(blocks)
            key=canonical(design)
            print(i,mode,key,sep='\t',file=out)
            counts[key]+=1
            if key in designs:
                continue
            name=f'design{len(designs):02}'
            F=base.frame(tiles)
            z=np.zeros(183,dtype=np.int64)
            z[14]=2
            for t in range(28):
                z[15+6*t:17+6*t]=1
            w=F.T@z
            assert set(w)=={-2,2}
            S=F.T@F-15*np.eye(496,dtype=np.int64)
            A=(S*w[:,None]*w[None,:]<0).astype(np.int64)
            base.verify(A,231,102,112)
            G=nauty_graph(A)
            (HERE/f'{name}.blocks').write_text(''.join(' '.join(map(str,b))+'\n' for b in blocks))
            (HERE/f'{name}.starter').write_text(''.join(' '.join(map(str,b))+'\n' for b in rows))
            (HERE/f'{name}.g6').write_bytes(graph6(A))
            designs[key]=(name,group_order(design),canonical(G),group_order(G))
    for k,(name,aut,h,g) in designs.items():
        print(name,'old design',k==original,'starters',counts[k],'design group',aut,'graph group',g,'graph certificate',h)
    print('distinct designs',len(designs),'total normalized starters',sum(counts.values()))


if __name__=='__main__':
    main(Path(sys.argv[1]))

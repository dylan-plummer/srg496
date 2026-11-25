"""Construct the graph on 976 points, its induced 496-vertex graph, and their descendants."""
from itertools import combinations,permutations
from pathlib import Path
import json
import networkx as nx
import numpy as np

HERE=Path(__file__).resolve().parent
START=HERE.parent/'starters'
W=np.array([[(-1)**((x&c).bit_count()%2) for x in range(8)] for c in range(8)],dtype=np.int16)
EYE=lambda n:np.eye(n,dtype=np.int16)


def tiles(design):
    rows=[list(map(int,l.split())) for l in (START/f'design{design:02}.starter').read_text().splitlines()]
    def move(s,i,j):
        return (56 if s==2 else 28*s+4*((i//4+j//4)%7)) + ((i^j)&3)
    return [{move(s,i,j):(a,b) for a,b,s,i in rows} for j in range(28)]


def profiles():
    """Choose a relative row permutation for each possible intersection dimension."""
    words=lambda order:{tuple(W[i,x] for i in order) for x in range(8)}
    first=words(range(1,8))
    ans={}
    for p in permutations(range(1,8)):
        h=len(first&words(p))
        ans.setdefault(h,p)
        if len(ans)==4:
            break
    assert set(ans)=={1,2,4,8}
    return {28+16*(k-1):v for k,v in ans.items()}


def make(d, record, colours):
    ts=tiles(d)
    frame=np.zeros((183,976),dtype=np.int16)
    for b,word in enumerate(record['matrix'][1:]):
        frame[b,:16]=[1-2*((word>>x)&1) for x in range(16)]
    rel=profiles()
    for p,colour in enumerate(colours):
        minus=slice(16+8*p,24+8*p)
        plus=slice(496+8*p,504+8*p)
        frame[p//4,minus]=1
        frame[p//4,plus]=-1
        flags=[(i,t[p]) for i,t in enumerate(ts) if p in t]
        cut=lambda f:(f[1][0]<2)!=(f[1][1]<2)
        allflags=sorted(flags,key=cut)
        assert len(allflags)==7 and sum(map(cut,flags))==(4 if p<56 else 0)
        extra=rel[colour] if isinstance(colour,int) else colour
        assert sorted(extra)==list(range(1,8))
        for c,(j,(a,b)) in enumerate(allflags,1):
            eps=((-1 if c==7 else 1)*(1 if a<2 else -1)) if cut((j,(a,b))) else 1
            u,v=15+6*j+a,15+6*j+b
            frame[u,minus]=eps*W[c]
            frame[v,minus]=-eps*W[c]
            frame[u,plus]=frame[v,plus]=W[extra[c-1]]
    assert np.array_equal(frame@frame.T,80*EYE(183))
    y=(frame.T@frame)-15*EYE(976)
    assert np.all(abs(y+EYE(976))==1)
    kernel=np.ones(976,dtype=np.int16)
    for i in (16,496):
        for p in range(60):
            if p%4<2:
                kernel[i+8*p:i+8*p+8]=-1
    assert not (frame@kernel).any() and not(frame[:,:496]@kernel[:496]).any()
    parent=y*kernel[:,None]*kernel[None,:]
    low=frame[:,:496]
    u=np.zeros(183,dtype=np.int16)
    u[14]=2
    for t in range(28):
        u[15+6*t:17+6*t]=1
    twice=u@low
    assert (abs(twice)==2).all()
    word=twice//2
    assert np.array_equal(y[:496,:496]@word,33*word)
    open_mat=y[:496,:496]*word[:,None]*word[None,:]
    def root(s):
        signs=s[0,1:]
        return s[1:,1:]*signs[:,None]*signs[None,:]<0
    return {'976_480':parent>0,'975_462':root(y),'496_240':parent[:496,:496]>0,
            '496_231':open_mat<0,'495_238':root(y[:496,:496])}


def export(d,h,colour,tag):
    m=make(d,h,colour)
    for name,a in m.items():
        g=nx.from_numpy_array(a)
        blob=nx.to_graph6_bytes(g,header=False)
        path=HERE/'graphs'/f'{tag}_{name}.g6'
        path.write_bytes(blob)
        print(path.name,len(g),min(dict(g.degree()).values()),max(dict(g.degree()).values()),flush=True)


if __name__=='__main__':
    (HERE/'graphs').mkdir(exist_ok=True)
    hh=json.loads((HERE/'hadamards.json').read_text())
    for d in range(3):
        for h in hh[1:]:
            export(d,h,[140]*60,f'd{d}_h{h["pointed"]}')
    systems={}
    for perm in permutations(range(1,8)):
        system=tuple((i,j,k) for i,j,k in combinations(range(7),3) if perm[i]^perm[j]^perm[k]==0)
        systems.setdefault(system,perm)
    assert len(systems)==30
    export(0,hh[4],list(systems.values())*2,'mixed_all')
    for k in (1,2):
        colours=[(44,76,140)[(p if k==1 else (p//4))%3] for p in range(60)]
        export(0,hh[1],colours,f'extension{k}')

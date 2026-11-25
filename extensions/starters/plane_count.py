"""Count independently chosen line structures and phases on the recovered fibres."""
from fractions import Fraction
from itertools import combinations,permutations
from math import factorial,prod
from classify import HERE,tiles_from
from recover import test


def plane_options(tiles,p):
    ts=[t for t,T in enumerate(tiles) if p in T]
    triples=set()
    for phi in permutations(range(1,8)):
        if p<56 and any((int(min(tiles[t][p])<2<=max(tiles[t][p])))!=(c>>2) for t,c in zip(ts,phi)):
            continue
        triples.add(tuple(sorted(sum(1<<ts[k] for k in tri) for tri in combinations(range(7),3)
                                 if phi[tri[0]]^phi[tri[1]]^phi[tri[2]]==0)))
    assert len(triples)==(6 if p<56 else 30)
    return triples


def count(path):
    rows=[list(map(int,l.split())) for l in path.read_text().splitlines()]
    tiles=tiles_from(rows)
    opts=[plane_options(tiles,p) for p in range(60)]
    choices=prod(map(len,opts))
    star=factorial(15)//20160
    # The intrinsic affine planes recover all cells and their direction data.
    plain=(choices*star+83)//84
    phases=2**(6*56+7*4)
    translations=16*8**60
    pivots={}
    for p in range(56):
        crossing=sum(1<<t for t,T in enumerate(tiles) if p in T and min(T[p])<2<=max(T[p]))
        assert crossing.bit_count()==4
        while crossing:
            pivot=crossing.bit_length()-1
            if pivot not in pivots:
                pivots[pivot]=crossing
                break
            crossing ^= pivots[pivot]
    tile_flips=2**(28-len(pivots))
    print(path.stem,'binary rank of crossing constraints',len(pivots))
    rich=(choices*star*phases+84*translations*tile_flips-1)//(84*translations*tile_flips)
    print(path.stem,'ordinary plane choices',choices,'star choices',star)
    print('lower bound from plane data',plain,'with free phases',rich)
    return plain,rich


if __name__=='__main__':
    counts=[count(p) for p in sorted(HERE.glob('design*.starter'))]
    print('total lower bounds',*(sum(v[i] for v in counts) for i in (0,1)))

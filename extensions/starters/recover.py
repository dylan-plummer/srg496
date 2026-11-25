"""Check intrinsic recovery of the fibres, independently of frame decorations."""
from collections import Counter
from itertools import combinations
from classify import HERE,tiles_from


def test(path):
    rows=[list(map(int,l.split())) for l in path.read_text().splitlines()]
    tiles=tiles_from(rows)
    flags=[{} for p in range(60)]
    for t,T in enumerate(tiles):
        for p,e in T.items():
            flags[p][t]=sum(1<<q for q,f in T.items() if q!=p and set(e)&set(f))
    assert all(len(s)==7 and sum(v.bit_count() for v in s.values())==56 for s in flags)
    # For a test fibre q, a point p has zero slope within its own part,
    # or contributes the character of the unique tile which joins p and q.
    slopes=[[sum(1<<t for t,neighbours in flags[q].items() if neighbours>>p&1)
              for q in range(60)] for p in range(60)]
    for size in range(1,5):
        survivors=[]
        for pts in combinations(range(60),size):
            dirs=[slopes[p] for p in pts]
            for col in zip(*dirs):
                what=0
                for val in col:
                    what ^= val
                if what.bit_count() in (1,2):
                    break
            else:
                survivors.append(pts)
        assert all(len({p//4 for p in pts})==1 and len(pts)>=3 for pts in survivors)
        print(path.stem,'odd ordinary fibres',size,'survivors',len(survivors),survivors[:6],flush=True)
    supports={}
    repeat=0
    group_unions=0
    for p,row in enumerate(flags):
        for fs in combinations(row.values(),4):
            m=fs[0]|fs[1]|fs[2]|fs[3]
            assert m.bit_count()==32
            if m in supports:
                repeat+=1
            supports[m]=p
            group_unions+=all(((m>>(4*g))&15) in (0,15) for g in range(15))
    assert not repeat and not group_unions
    print(path.stem,'two-pair support collisions',repeat,'part unions',group_unions,flush=True)

if __name__=='__main__':
    for path in sorted(HERE.glob('design*.starter')):
        test(path)

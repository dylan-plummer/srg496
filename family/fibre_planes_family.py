"""Run the zero-sum 4-set fibre test on random family members."""
import sys, random
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'audit'))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from independent_rebuild import tiles_from, read_starter
from family import member
from fibre_planes import zero_sum_quads, fibre

N = int(sys.argv[1]) if len(sys.argv) > 1 else 5
rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 7)
tiles = tiles_from(read_starter())
for i in range(N):
    A = member(tiles, rng)
    Q = zero_sum_quads(A)
    inside = sum(len({fibre(v) for v in q}) == 1 for q in Q)
    print(f'member {i}: {len(Q)} zero-sum 4-sets, {inside} inside one fibre', flush=True)

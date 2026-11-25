# Three supporting designs and a lower bound for the 496-vertex row

24 November 2026

There are at least

\[
  \frac{6^{56}30^4}{84}\frac{15!}{20160}\,2^{152}(2^{22}+2\cdot2^{19})
  > 7.06\times 10^{107}
\]

pairwise nonisomorphic strongly regular graphs with parameters
\((496,231,102,112)\). In fact this bounds the number of distinct regular
two-graphs in the construction. The same lower bound therefore holds separately
for the degree-264 complements, for \((495,238,109,119)\), and for
\((496,240,120,112)\). The last two parameter sets already have known examples.
The bound counts constructed objects up to isomorphism. It makes no assertion
that all its members have been individually exported or compared to the
literature.

## The three starters

The files in [starters](starters/) give three decompositions of the complete
15-partite graph on 60 points, four per part, into 28 triangular graphs
\(T(6)\). The point and group numbering and the frame construction are those
of the [original proof](../seed/PROOF.txt).
Each starter file contains 15 rows giving the two edge symbols, the sheet,
and the position. Translating the 15 entries produces every tile.
For each tile the six symbol stars give six blocks of size five. Adjoin
point 60 to each of the 15 parts for the other 15 blocks. Each of the three
block files specifies a Steiner \(2\!-(61,5,1)\) design.

The normalized starter search returned 216 certificates. Canonicalizing the
uncolored point-block incidence graphs gave three classes of 72. Their
automorphism groups all have order 84. Points and blocks have different
incidence degrees, 15 and 5. One design is the earlier seed. The other
two give designs and degree-231 graphs different from the earlier seed.
These design classes are relative to this calculation. See the
[starter census](starters/results.txt).

For every starter the raw frame has 16 special columns and eight columns
at each of 60 ordinary points. Its squared column norm is 15 and its
frame bound is 48. On each long-sheet point the four crossing tile
flags are the complement of a line in the Fano plane of nonzero Walsh
characters. Choose the signs on these four roots to have product −1
in the orientation toward the marked side of the tile. Quiet flags
have arbitrary signs. The special and short-sheet terms and the
four-character identity in the original proof give a sign eigenvector
of the raw Seidel matrix \(S\) with eigenvalue 33. The resulting regular
switch has degree 231 for every decoration used below.

## Recovering the fibers under arbitrary Seidel switching

Put \(M=S-I\), a sign matrix with diagonal −1. Call a quadruple of
distinct vertices admissible if the product of their four rows of \(M\)
is a constant vector. This property survives every simultaneous row and
column switch and every relabeling. Take connected components of the
hypergraph of admissible quadruples.

Inside each ordinary cell of eight vertices and the special cell of
sixteen, the admissible quadruples are exactly the affine planes on its
binary coordinate space. The row product evaluated on other cells tests
all nonzero Walsh characters, and within the cell every entry of \(M\)
is −1. Thus each cell is connected in the hypergraph.

The following finite checks depend only on the small unsigned starters.
They hold for every admissible choice of labels and phases.

1. For a set \(O\) of ordinary cells occurring an odd number of times
   in a quadruple, examine its row product on an ordinary cell \(q\).
   A member of \(O\) in the part of \(q\) contributes zero Walsh slope.
   Every other member contributes the character of the unique tile
   connecting it to \(q\). Take these tile flags with multiplicities
   modulo two. One or two surviving distinct flags give nonzero slope
   for every labeling. The special cell has zero slope on \(q\).
   Checking all subsets of at most four ordinary cells leaves no one-
   or two-element sets. Every surviving triple or quadruple is contained
   in one part. A surviving triple must be accompanied by one special
   vertex. On the special cell it then has the nonzero character of that
   part, so it fails there. A surviving quadruple uses one vertex in
   each of the four cells of one part. Its product on the special cell is
   the fourth power of that part character, hence +1. On each of its own
   four cells exactly one row contributes −1 and the others +1. It
   therefore fails as well.
2. When \(O\) is empty, a possible mixed quadruple uses two pairs in
   distinct cells. The difference of two rows of an ordinary cell,
   in binary notation, is zero on the special cell and is supported on
   four of its seven sets of eight neighboring ordinary cells. These
   sets are disjoint. Across all 60 points and all choices of four flags,
   the resulting 2,100 supports are distinct. None is a union of the
   four-point parts. The difference of a pair of special rows is
   supported on a union of parts and is zero on the special cell.
   Consequently mixed two-pair quadruples always fail. The common value
   on the special cell fixes the otherwise possible complement ambiguity.

The enumerations take only 60 points as input.
[recover.py](starters/recover.py) checks these claims with assertions.
Its [output](starters/recovery.out) gives all survivor counts. The admissible
quadruples are therefore *exactly* the 140 affine planes of the special
cell and the 14 affine planes in each of the 60 ordinary cells, 980 in
all. Their hypergraph components recover all 61 cells, with the special
one distinguished by its size. This holds for every permitted decoration.

On the special cell the projective row patterns of all ordinary vertices
are the 15 distinct part characters. Here projective means up to a
constant sign. Switching multiplies every pattern by a common sign
function on the cell. Equality of the patterns recovers the 15 parts.

On an ordinary cell \(p\), the special vertices provide a common projective
row direction. Dividing all external directions by it removes the switching
on \(p\). Each other part then exposes the Walsh character of the flag that
connects its point to \(p\). Form vertices \((p,c)\) for the seven flags at
each point. For every pair of base points in different parts, join the flags
they assign to one another. The connected components are exactly the 28
copies of \(T(6)\). The 5-cliques of each component recover its six symbol
stars. This reconstructs the entire supporting Steiner design. An independent
[blind audit](starters/blind_audit.py) starts with each exported graph, switches
and relabels it, and recovers its original design by this procedure.

## Counting inequivalent decorations

On a labeled seven-element flag set there are \(7!/168=30\) Fano plane
structures. Requiring the three quiet flags to be a line leaves
\(30\cdot 7/35=6\) choices. This gives \(6^{56}30^4\) independent choices
on the ordinary points. On the 15 parts one may choose any projective
3-space structure over \(\mathbb F_2\), giving \(15!/20160=64,864,800\)
choices. Fix one representative Walsh labeling for each structure.

All those structures are intrinsic to the two-graph together with the
recovered base. On an ordinary cell the triple product of normalized
projective character directions recognizes the Fano lines. On the
special cell the 15 directions are a translate of the 15 nonzero
characters, in the multiplicative group of sign functions up to
constant sign. Their affine span has size 16. The unique missing
direction identifies the origin, so the PG(3,2) incidence is recovered
by multiplication after translating the origin back.

For fixed labels, phases have \(2^{6\cdot56+7\cdot4}=2^{364}\) choices.
Each long point has three free quiet signs and three free crossing
signs. Each short point has seven free signs.

Fix a source two-graph on one of the three designs. A target isomorphic
to it induces one of at most 84 design permutations, since the design is
intrinsic. Each such permutation determines the target Fano and projective
structures and consequently their chosen Walsh labelings. The affine
planes describe the affine structure in each cell. Once the directions
carried by an isomorphism have been prescribed, the linear part of the
isomorphism on that cell is determined. The remaining freedoms are
translations, at most \(16\cdot8^{60}=2^{184}\) choices for the entire
vertex permutation. Equivalently, the ratios of the character values at
two vertices prescribe their binary difference.

For any fixed labeling and fixed vertex permutation, at most one binary
choice per tile can identify target phase assignments in a labeled switching class.
Indeed two raw target matrices with the same labels agree within the
special cell and between that cell and every ordinary vertex. A switching
relating them must first be constant on the special cell and then on
every vertex. Thus the two raw matrices agree entry by entry. On an
edge of a tile this fixes the product of the two flag phases. Since
\(T(6)\) is connected, all the phases in each tile agree up to one
common reversal. At each of the 56 long points, the product of the four
crossing phases must remain −1. The subset of reversed tiles consequently
satisfies 56 homogeneous binary parity equations. The 56-by-28
crossing incidence matrix has rank 22 on the first design and 19 on
each of the other two. Hence at most \(2^6,2^9,2^9\) tile reversals,
respectively, preserve the admissibility of target phases. These ranks
are computed directly from the fifteen-entry starters in plane_count.py.
Even allowing all tile reversals without this refinement would already
give more than \(4.03\times10^{101}\) distinct two-graphs.

Dividing the number of choices by these multiplicities and adding the
three disjoint design classes gives the displayed lower bound.
[plane_count.py](starters/plane_count.py) independently enumerates the
local Fano structures and evaluates the integer expression.
Every two-graph yields a degree-231 regular member. Isomorphic graphs
necessarily yield isomorphic two-graphs.

## The other two parameter sets

Isolate one vertex by Seidel switching and delete it. Its row in the
switched Seidel matrix is positive off the diagonal. The identity
\(S^2=18S+495I\) gives row sum 18 on the remaining 495 vertices.
Their adjacency graph has parameters \((495,238,109,119)\). The
descendant determines its parent two-graph by adjoining an isolated
vertex, so different two-graphs yield different descendants.

For a degree-240 member choose the positive constant sign word on
the special cell. Within each part choose two of its four ordinary
cells to carry the constant sign −1, and the other two +1. This word
\(w\) satisfies \(Fw=0\). The group-coordinate sums cancel, every
nonconstant ordinary Walsh character sums to zero, and the special
columns sum to zero. Switching by \(w\) gives Seidel row sum −15.
The positive off-diagonal entries give the degree-240 graph with
parameters \((496,240,120,112)\). Its complement has degree 255.
Consequently the lower bound from the three designs transfers to this
parameter set as well.

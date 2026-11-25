# Transfers to other rows

## Tight subframes at 196

No member of the 496-column construction on these starters contains a
196-column subframe whose Seidel eigenvalues in the raw orientation are
13 and −15. This is the orientation that would produce the degree-90
row by regular switching and complementation. More generally the
argument uses only the group-coordinate rows and the Walsh simplex
at the special cell.

Suppose a subset of 196 columns is tight with frame bound 28, as those
eigenvalues require. Let \(s\) be the number of selected special columns,
\(n_g\) the number of selected ordinary columns in part \(g\), and
\(t_p\) the number at its ordinary point \(p\). For a selected column
\(f_i\), the partial frame operator satisfies
\(\sum_{j\in X} f_j\langle f_j,f_i\rangle=28f_i\).
At the group coordinate of a selected ordinary column, special columns
contribute \(s\), the other cells in its part contribute \(n_g-t_p\),
and its own cell contributes \(16-t_p\). Thus
\[
              n_g+s=12+2t_p\quad(t_p>0).
\]
Every occupied cell of a part has the same occupancy. If \(a_g\) of its
four cells are occupied then \((a_g-2)t_p=12-s\).

If \(s=0\), each occupied part must use four cells of size six. The total
196 is indivisible by 24. If \(s>0\), apply the same frame identity to
a selected special vertex \(y\) and group coordinate \(g\). With \(R\)
the selected set in the binary 4-space one obtains
\[
  \sum_{z\in R}\chi_g(z)=(n_g-12)\chi_g(y).
\]
Every character varying on \(R\) has zero sum there and \(n_g=12\).
Every character constant there gives \(n_g=12+s\). The Fourier support
of the indicator of \(R\) therefore lies in the annihilator of the
span of its differences. Fourier inversion makes the indicator
constant on cosets of that span. As \(R\) lies in one such coset, it is
the whole coset. Its size can only be 1, 2, 4, 8, or 16.

| Special size | Forced occupancy contradiction |
|---|---|
| 1 | Every part would have 13 points drawn equally from at most four cells of capacity eight. |
| 2 or 4 | A varying character gives \(n_g=12\) and cell size \(s/2\), requiring more than four cells. |
| 8 | The nontrivial constant character gives \(n_g=20\) and cell size eight. |
| 16 | Every part gives \(n_g=12\) and cell size eight. |

The opposite Seidel orientation for the 196 row would have the two raw
eigenvalues 15 and −13. A 196-column subset of a rank-155 frame
necessarily retains the eigenvalue −15, so this orientation is also
impossible. The finite order-seven model in transfer196.py corroborated
the first orientation on each of the nine symmetric examples.

There is also a size obstruction to tight subframes of orders 376 and
406. A proper tight subframe with frame bound below 48 has complement
spanning the full 155-dimensional space, because its frame operator is
\(48I-LP\) for an orthogonal projector \(P\) and \(L<48\). The two
complements would have only 120 and 90 columns.

## Two other attempts

At 148, a mixed template uses one A6 tile, three E8 tiles and seven
group coordinates. Its sign word for degree 63 would force a square
congruent to 3 modulo 7 in the discriminant form of A6. This eliminates
all placements and Walsh labelings of these tiles. See the
discriminant argument (work directory `mixed148_20261124`, not included here).

At 288 the proposed base has six parts of twelve, four columns per base
point, and three spanning root tiles. A tile takes both root lines
\(e_i+e_j\) and \(e_i-e_j\) over every edge of a simple 6-regular graph
on twelve symbols. Its 72 lines have frame operator 12 and every line
has twenty nonorthogonal neighbors. A decomposition of the 6-partite
base would therefore produce a frame of rank 42, norm seven and bound
48. Giving every group coordinate the value one supplies a sign word
for degree 123. Balancing ordinary cell signs in each part supplies a
kernel word giving the degree-140 complement. Isolating a vertex would
give the 287-vertex descendant.

I tested regular actions on the 72 ordinary points which fix the three
tiles and act through permutations of their twelve coordinate symbols.
For a connected or disconnected edge-transitive spanning graph, the
symbol action has either one orbit, or two endpoint orbits of size six
each. The latter case gives the complete bipartite graph. The former
uses the index-12 subgroup for a symbol and an unordered pair with
orbit size 36. The acting group is regular on the base points, so
neighborhoods of the identity determine all edges. The size-12 part
through the identity must be a subgroup. The three tile neighborhoods
must be disjoint 20-subsets partitioning its complement.
`regular_tiles.py` and
`groups.g` (work directory `root_tiles_288_20261124`, not included here) checked the 50 groups
of order 72 and found no such partition. Actions permuting the tiles
and more general placements remain open.

## Involutions in the one-starter approach at 616

Take a single T(7) starter developed by a group of order 25 on three
regular sheets of 25 ordinary points each. Each sheet is partitioned
into the five cosets of a subgroup of order five. A point appears in
seven tiles, so the starter contains seven points of each sheet. An
involution normalizing the translation group would permute the 25 tiles
and hence fix a starter tile. Write \(M\) for the induced automorphism
of the group of order 25. Origins can be chosen on invariant sheets at
fixed points and compatibly on exchanged sheets. The involution then
acts linearly by \(M\) between sheets. This follows also by halving
the translation offset in the odd-order group.

A same-sheet translation class is a difference outside the part subgroup,
up to sign. A cross-sheet class is an oriented difference in the full
group. The fixed starter would contain exactly one edge per class, so
the number of fixed classes must equal the number of fixed edges of
its T(7). On the seven symbols an involution with one, two, or three
transpositions fixes respectively 35, 9, or 3 edges of T(7). The
identity fixes 105. These counts follow by taking unordered pairs of
intersecting two-subsets of the seven symbols.

The possibilities for \(M\) are the identity, negation, and, only on
the elementary abelian group, a reflection with two eigenlines. On a
fixed sheet these fix respectively 10, 10, and 2 same-sheet difference
classes. In the reflection case the invariant part subgroup must be
one of the two eigenlines.

| Action on sheets | Identity on translations | Negation | Reflection |
|---|---:|---:|---:|
| All three fixed | 105 | 33 | 21 |
| Two exchanged | 11 | 35 | 7 |

For three fixed sheets the cross classes contribute three copies of
\(|\ker(M-I)|\). For exchanged sheets the only additional fixed classes
are between the two exchanged sheets and contribute \(|\ker(M+I)|\).
The 105 case is the identity action on all base points. In the only other
numerical match, negation with two sheets exchanged, all 25 cross
classes between those sheets are fixed. A fixed edge between exchanged
sheets must exchange its endpoints. The starter has only seven points
on either sheet and can supply at most seven such edges. This rules out
every nontrivial involution normalizing these one-starter developments.
The underlying 616 construction with odd symmetry remains open.

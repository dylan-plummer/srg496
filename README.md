# Strongly regular graphs with parameters (496, 231, 102, 112)

Graphs, code and a short note (`note/srg496.pdf`) for the results below (Dylan Plummer, 24 November 2026).

## Results

| parameters | status in Brouwer's table (checked 24 Nov 2026) | here |
|---|---|---|
| srg(496, 231, 102, 112) and complement (496, 264, 144, 136) | open (`?`, "2-graph?") | more than 3.43 · 10^110 pairwise non-isomorphic graphs (lower bound, three designs, five Hadamard matrices) |
| srg(496, 240, 120, 112) and complement (496, 255, 126, 136) | known (Wallis; NO⁺(10,2)) | the same lower bound; exactly 554,146,599,744 classes in one explicit subfamily |
| srg(495, 238, 109, 119) | known (O⁻(10,2), Mathon) | the same lower bound (descendants) |
| srg(976, 480, 248, 224) | known (Wallis) | more than 5.04 · 10^86 non-isomorphic graphs containing one fixed labeled srg(496,240,120,112) as an induced subgraph |

- **Tiling.**
  - One starter under Z7 × 𝔽₂² decomposes K_{15×4} into 28 copies of T(6). Two sheets have 28 points, one has 4,
    and the parts have size 4.
  - Adjoining one point gives a Steiner system S(2,5,61) (`data/steiner_2_61_5.txt`).
- **Frame.**
  - Coordinates: 15 group coordinates and one ℝ⁶ per tile.
  - The added point gets 16 Walsh columns, and every other point gets 8, one for each x ∈ 𝔽₂³. Each column holds
    a 1 in its group coordinate and χ_c(x)·(e_a − e_b) in each of its 7 tiles, with the 7 nonzero characters
    assigned to the 7 tiles.
  - This is a real equiangular tight frame: 496 vectors of rank 155, Seidel spectrum 33^155 (−15)^341.
- **Sign vector.** Mark symbols 0, 1 of every K6.
  - Every long-sheet point crosses the cut in 4 tiles. Give those tiles the 4 characters outside a hyperplane
    of 𝔽₂³, with sign product −1.
  - Then w = Fᵀz, with z = ½ on the marked coordinates and 1 on the short sheet's group coordinate, is a ±1
    vector with Sw = 33w.
  - Switching by w gives the regular graph of degree 231, with spectrum 231^1, 7^341, (−17)^154.
- **Exported graph.** |Aut| = 4 (C2 × C2), checked with both nauty and bliss.
- **Steiner half-frame.** The frame is exactly half of the Steiner ETF of the S(2,5,61): 976 vectors in ℝ¹⁸³,
  Seidel spectrum 65^183 (−15)^793 (`audit/steiner_subframe.py`). So the two-graph is a regular sub-two-graph of
  a Steiner two-graph.
  - This does not contradict the corollary in our 442 paper. That corollary is about the Steiner two-graphs
    themselves: on 496 points from an S(2,3,31), and on 976 points from an S(2,5,61).
- **Two-graph.** Ω is not the symplectic regular two-graph on 496 points, nor its complement (`audit/twograph.py`).
  - The symplectic one is 2-transitive, with descendant |Aut| = 5.003·10¹³ = |Sp(10,2)|/496.
  - Ω's descendants have |Aut| = 4 at the special columns and 1 at ordinary ones, so Aut(Ω) is not even
    transitive.
  - By the corollary above, Ω is also not the Steiner two-graph of any S(2,3,31).
  - Its descendants are srg(495,238,109,119) (a known row) with at most 4 automorphisms. We make no novelty
    claim for Ω itself.
- **Family.** Several choices are free:
  - at each long-sheet point: the hyperplane, both orders of characters, and the signs (64,512 ways);
  - at each short-sheet point: the labels;
  - at the added point: the parts-to-characters bijection.

  Fifty random members were pairwise non-isomorphic, each with a trivial automorphism group
  (`family/family_sample.out`).
- **Three designs** (`extensions/`, companion work). A search for starters that also have an order-3
  symmetry found 216 starters in three design classes. Each S(2,5,61) has |Aut| = 84, and the three are
  pairwise non-isomorphic. One is the original design.
- **Recovery from the two-graph.** For every admissible decoration of the three designs and any Seidel
  switching, the 4-sets whose rows of S − I have constant product are exactly the 980 affine planes inside
  cells (finite check on the starters). From them the two-graph determines the cells, parts, tiles, design,
  the Fano structure at each point and the PG(3,2) structure on the parts.
- **Lower bound** (note, Theorem 7). There are more than
  6^56 · 30^4 / 84 · 15!/20160 · 2^152 · (2^22 + 2·2^19) > 7.06 · 10^107 pairwise non-isomorphic regular two-graphs.
  - Each one gives an srg(496,231,102,112) and, through a kernel word constant on cells, an
    srg(496,240,120,112). Its descendants are srg(495,238,109,119). So the bound holds for all three rows.
  - The 2^152 · 2^22 (or 2^19) factor comes from the signs, divided by the 2^184 cell translations and by the
    tile reversals allowed by the crossing matrix (GF(2) rank 22 on the first design, 19 on the others).
  - Without that last refinement the bound is still 4.03 · 10^101.
  - Every member's automorphism group fixes the 16-point special cell, so no member's two-graph is transitive.
- **Other Hadamard matrices (25 Nov).**
  - The 16 columns at the added point can come from any normalized Hadamard matrix of order 16. There are five up
    to equivalences fixing the all-ones row.
  - Their special-cell four-set counts are 140, 76, 44, 28 and 28. Their column groups have orders 322560, 9216,
    1536, 2688 and 2688.
  - The bound grows by the factor 1 + 35 + 210 + 120 + 120 = 486, to more than 3.43 · 10^110 (note,
    Proposition 5.3).
- **Steiner completion (25 Nov).**
  - At each point, add a second group of 8 columns, with −1 on the part row and χ(x)(e_a + e_b) in the tiles.
    This gives the full Steiner frame on 976 points, with Seidel eigenvalues 65 and −15.
  - The kernel word gives srg(976,480,248,224), whose first 496 vertices induce the srg(496,240,120,112).
    Descendants are srg(975,462,205,231). Both rows are known.
  - The 30^60 choices of the second Fano planes give at least ⌈30^60/84⌉ > 5.04 · 10^86 inequivalent parents of
    one fixed labeled half (note, Proposition 8.1).
  - Data and checks are in `extensions/hadamard/`.
- **Symmetric members.** On the first design, nine two-graphs have automorphism group of order 84.
  - They give 56 pairwise non-isomorphic srg(496,231,102,112) by cuts, with |Aut| in {2, 12, 14, 84}.
  - They have 126 descendant classes, with |Aut| in {1, 3, 12, 21, 84}.
  - Their cell-constant kernel words give exactly 554,146,599,744 classes of srg(496,240,120,112), by Burnside.
- **Template.**
  - The same frame works for any T(m)-decomposition of K_{P×(m−2)}, given Hadamard matrices of orders P+1 and
    (P+1)/2 (`template/scan_template.py`).
  - Up to 1300 vertices it reaches known rows: 10, 16, 36, 64, 136, 256, 528, 820 and 1024.
  - (148 would need an S(2,6,36), the nonexistent affine plane of order 6.)
  - Besides 496 the only open targets are on 615 and 616 vertices. They need an S(2,6,76) whose blocks off one
    point split into 25 copies of T(7).
  - That tiling alone would give srg(615,294,133,147) as a descendant. A sign vector would give
    (616,287,126,140) or (616,315,154,168).
  - Neither published S(2,6,76) admits such a tiling.

## Files

- `data/`: `starter.tsv`, `steiner_2_61_5.txt`, `brouwer_table_snapshot.csv`.
- `graphs/`: `srg_496_231_102_112.g6`, `srg_496_264_144_136.g6`, `SHA256SUMS`. Vertices: 16 special columns,
  then 8 columns for each point 0..59.
- `seed/`: the original material.
  - `PROOF.txt`: the proof.
  - `construct.py`: exact integer constructor.
  - `audit.py`: standard-library all-pairs audit plus an edge-by-edge rebuild from the starter.
  - `symmetry.py`: nauty and bliss.
  - `search_starter.py`: regenerates the starter.
- `audit/`:
  - `independent_rebuild.py`: rebuild from the written description only, checked for byte equality with the
    exported graphs.
  - `steiner_subframe.py`: the Steiner half-frame.
  - `twograph.py`: the two-graph comparison.
- `family/`:
  - `family.py`: random members, srg check, nauty certificate and |Aut|.
  - `fibre_planes.py`, `fibre_planes_family.py`: fiber recovery.
- `template/`: `scan_template.py` and its output.
- `extensions/`: the companion analysis (snapshot of 24 Nov 2026) plus our audits.
  - `PROOF.md`: the lower bound. `CATALOGUE.md`: the order-84 seeds and the reproduction commands.
    `TRANSFERS.md`: attempts at the 196, 148, 288 and 616 rows (not audited in depth).
  - `starters/`: the three starters with their block files, the starter census (`enumerate.cpp`,
    `classify.py`), the recovery checks (`recover.py`, `blind_audit.py`) and `plane_count.py`.
  - `graphs/`: the 56 degree-231 graphs and 126 descendants (`all_*.g6.gz`), and individual seed graphs
    for the rows 231, 240 and 495. `SHA256SUMS` covers these files and the starters.
  - `family.py` and the other companion scripts: the seed constructor, the census, the switch classes and
    the Burnside count.
  - `designs_audit.py` (ours) checks the designs, the GF(2) ranks, random members of all three designs in all
    three rows, recovery under random switching, the Fano count and the bound.
  - `hadamard/`: the five Hadamard matrices, the 976-point completion, 75 example graphs, the companion audit and
    recovery checks, and our `independent_hadamard.py`.
  - `seeds_audit.py` (ours) checks, with nauty, that the raw group and the two-graph group have order 84, the 126
    descendant classes, every exported graph, and the Burnside count recomputed from the nauty group.
- `note/`: `srg496.tex` and `srg496.pdf`.

## Verification

    ./verify_all.sh            # about 2 minutes; PYTHON=... to choose the interpreter

Requirements: Python ≥ 3.10 with numpy, networkx, python-igraph and pynauty. The `seed/audit.py` check
uses only the standard library.

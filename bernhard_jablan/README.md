# The knot 13n1587 (Section 5)

For a nontrivial knot $K$, $N_{\min}(K)$ is the set of knots obtained by changing one crossing in a minimal diagram of
$K$, and $u^s_{\mathrm{BJ}}(K)=1+\min\{u(K'):K'\in N_{\min}(K)\}$. The paper shows $u(13n1587)=2$ and
$u^s_{\mathrm{BJ}}(13n1587)=3$.

## Minimal diagrams (`code/`, records in `results/`)

The minimal diagrams are enumerated following Brittenham–Hermiller, Prop. 2.3. A minimal diagram of $K$ arises,
up to flypes, from one reduced alternating reference diagram per alternating knot with the same crossing number,
by crossing changes.

| Script | Output | Content |
| --- | --- | --- |
| `enum_diagrams.py TARGET` | `diagrams_TARGET.json` | all diagrams obtained from the reference diagrams by crossing changes, filtered by the determinant and identified with SnapPy |
| `neighbors.py diagrams_TARGET.json` | `bjset_TARGET.json` | the knots obtained by one crossing change in each minimal diagram |
| `verify_completeness.py TARGET OUT.json` | `completeness_TARGET.json` | completeness with exact invariants (details below) |
| `cover_invariant.py TARGET OUT.json DMAX` | `cover_separation_TARGET.json` | separates candidates with the same Alexander polynomial and knot Floer homology, by the homology of all connected covers of degree at most DMAX |

`verify_completeness.py` works as follows:
- It computes the determinant of every candidate exactly, by fraction-free elimination.
- It compares exact invariants of every candidate with those of the target: the crossing number after
  simplification, the Alexander polynomial, and knot Floer homology up to mirror.
- It stores every compared value.

For 13n1587 the counts are:
- 4,878 reference diagrams;
- 13,346 diagrams with determinant 75;
- 11,903 simplify to fewer than 13 crossings;
- 1,379 have a different Alexander polynomial;
- the remaining 64 all have the Alexander polynomial and the knot Floer homology of 13n1587, so knot Floer
  homology excludes none of them; `cover_invariant.py` separates 57 of them by covers of degree at most 4
  (computed up to degree 5);
- the remaining 7 are minimal diagrams of 13n1587.

The 7 minimal diagrams give 91 crossing changes and 8 knot types. For 12n563, 1,288 reference diagrams give 18
minimal diagrams. None of their 216 crossing changes gives the unknot (`bjset_12n_563.json`).

## Lower bounds for the neighbours

| Script | Output | Content |
| --- | --- | --- |
| `neighbor_lower.py bjset_13n_1587.json` | `neighbor_lower_13n_1587.json` | $u\ge2$ from the diagram: $\lvert\sigma\rvert\ge4$, two generators of $H_1(\Sigma)$ (Wendt), or the rank-one test |
| `composite_neighbors.py` | `composite_neighbors_13n_1587.json` | the composite neighbours: signature and summands |
| `neighbor_rigor.py` | `neighbor_rigor_13n_1587.json` | cover homology of each neighbour against its named knot; for 10₃₈, a reduced alternating diagram in which no crossing change gives the unknot (McCoy's theorem) |
| `cert_identify.py`, `cert_identify_deep.py`, `cert_identify_deep_generic.py` | `cert_identify_13n_1587.json` | identification of 11n66, 11n99 and 12n563 by combinatorial isomorphisms of triangulations of the complements, found independently with SnapPy and with Regina |
| `hfk_crosscheck.py` | `hfk_crosscheck_13n_1587.json` | knot Floer homology of each neighbour against its named knot |
| `u1_independent.py NAME=PATH ...` | `u1_independent_*.txt` | independent Ni–Wu test of $u=1$ on given correction terms. The Greene-model correction terms of 11n66 and 11n99 are in `greene_u1_11n_*.json` |

`owens_independent.py` is an independent implementation of Owens's obstruction, used by `../plumbing/`. `bf2.py` is
the same file as `../rank_two/bf2.py`.

## Upper bound (`upper_bound/`)

`check.py` writes `check.json` with three parts:
- **(a)** `hits.jsonl` contains a 15-crossing diagram of 13n1587 in which changing two crossings gives the unknot.
  - The diagram is identified with KnotInfo's 13n1587 by isomorphic triangulations of the complements.
  - The changed diagram simplifies to 0 crossings, has knot Floer homology of total rank one, and has
    $\pi_1\cong\mathbb Z$.
  - The same search also found 8 diagrams of 13n1587 in which one crossing change gives $10_{113}$; they are
    identified in the same way.
- **(b)** The Casson–Walker identity holds for the Greene-model correction terms of 11n66, 11n99 and 12n563.
- **(c)** The plumbing route for the same knots found no negative-definite plumbing.

The other upper-bound argument of the paper uses the braid representatives of 13n1587 and $10_{113}$ of
Brittenham–Hermiller, which differ in one crossing. They are verified by `upper_bounds/verify_braid_certificate.py`
in the [arXiv:2609.09861 repository](https://github.com/Pinocchio315/unknotting-tait-graphs), with the same kind of
identification.

## Commands

Run from `code/`, writing to a new folder so that the deposited records stay unchanged.

```sh
OUT=$(mktemp -d)
python -B verify_completeness.py 13n_1587 $OUT/completeness_13n_1587.json
python -B cover_invariant.py 13n_1587 $OUT/completeness_13n_1587.json 5
python -B verify_completeness.py 12n_563 $OUT/completeness_12n_563.json
python -B enum_diagrams.py 13n_1587 && python -B neighbors.py diagrams_13n_1587.json
python -B neighbor_lower.py ../results/bjset_13n_1587.json
python -B composite_neighbors.py ../results/bjset_13n_1587.json $OUT/composite_neighbors_13n_1587.json
python -B neighbor_rigor.py ../results/bjset_13n_1587.json $OUT/neighbor_rigor_13n_1587.json
python -B cert_identify.py ../results/bjset_13n_1587.json 13n_1587 $OUT/cert_identify_13n_1587.json
cd ../upper_bound && python -B check.py
```

`verify_completeness.py` evaluates about twenty million diagrams. `neighbor_lower.py` takes seconds, and `check.py`
takes about 3.5 minutes.

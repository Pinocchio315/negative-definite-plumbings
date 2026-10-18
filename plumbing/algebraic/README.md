# Algebraic diagrams (Lemma 2.7, Sections 3.1 and 3.4)

For the 32 knots of `targets.txt`, Regina recognises neither the double branched cover nor the pieces cut along an
incompressible torus, and no piece is certified hyperbolic (`../jsj.json`). These knots have algebraic diagrams, and
Lemma 2.7 reads a plumbing tree of $\pm\Sigma(K)$ from such a diagram. The first version of the paper left these knots
undecided; the revised version places them in class (2).

Run the scripts from this folder. Each script writes the JSON file of the same name next to it.

| Script | Output | Content |
| --- | --- | --- |
| `basic_polyhedron.py` | | merging two vertices that bound a bigon, starting from a PD code (edges counterclockwise at each crossing). A diagram is algebraic if merging reaches a single vertex. `validate` compares this with the Conway notations in KnotInfo (at most 12 crossings) |
| `algebraic_plumbing.py` | | the tree of Lemma 2.7, by the rules (a)–(c) of Section 2.2, and its reduction to the Neumann normal form by Neumann moves |
| `alg_diagrams.py` | `alg_diagrams.json` | algebraic diagrams of the 32 knots, found by random Reidemeister moves from the KnotInfo diagrams. The search is seeded, but spherogram's moves are not exactly reproducible, so a rerun may find other diagrams; these give the same normal forms |
| `u_plumbing.py` | `u_plumbing.json` | for both orientations: the normal form, $\lvert\det Q\rvert$, the cokernel, the exact inertia and the nodes. Also the orientation from the Casson–Walker identity and the linking form, which the conclusions do not need |
| `crosscheck_u.py` | `crosscheck_u.json` | a second algebraic diagram of each knot (another seed) gives the same two normal forms. SnapPy finds both diagrams isometric to the KnotInfo diagram, and NumPy eigenvalues agree with the exact inertia |
| `validate_plumbing.py` | `validate_plumbing.json` | the construction on the 3,848 KnotInfo diagrams with at most 13 crossings that are algebraic: determinants, and agreement with the trees identified by triangulation in `../neumann/identify_trees.json` |
| `cw_general.py` | `cw_general.json` | the Casson–Walker formula with the signature term, $L=(\lvert\det Q\rvert/12)\bigl(\operatorname{tr}Q-3\sigma(Q)+\sum_v(2-\delta_v)(Q^{-1})_{vv}\bigr)$, which is the formula of Section 3.4 for negative-definite trees. It checks invariance under the Neumann moves used and finds the orientation with $L=R(K)$ |

## Results

- **The KnotInfo diagrams.** For the 32 knots, merging stops at six vertices. The algebraic diagrams found have 14 to
  16 crossings.
- **The 32 knots.** For each knot, the normal forms of $\Sigma(K)$ and $-\Sigma(K)$ satisfy the following.
  - Both have four nodes, so $\Sigma(K)$ has four Seifert pieces.
  - Both have positive eigenvalues: $p^+\in\{1,2,3\}$.
  - The cokernel of $Q$ is cyclic of order $\det K$.
  - A second algebraic diagram gives the same normal forms.

  Hence all 32 knots are in class (2). `../../paper/make_tables.py` asserts these conditions.
- **The construction on the KnotInfo diagrams.** On the 3,848 algebraic KnotInfo diagrams:
  - every determinant agrees;
  - for the 165 knots of class (2) among them, the two normal forms are the two trees identified by triangulation;
  - for the 151 knots of class (3), the negative-definite tree is one of the two normal forms.

  These 316 knots include 33 of the 43 knots of Tables 3 and 4.
- **Orientation.** For the KnotInfo diagrams, the sign $\varepsilon=-1$ of the construction gives $\Sigma(K)$. The
  Casson–Walker identity holds for this sign on 3,813 diagrams, and for both signs on the 35 with $L=R=0$.

## Rerun

On 19 October 2026 the five scripts were rerun in this layout.
- The search found other diagrams for 30 of the 32 knots, again with 14 to 16 crossings.
- The normal forms, their inertia, the Casson–Walker signs and all checks of `crosscheck_u.py` were the same.
- `validate_plumbing.json` and `cw_general.json` were reproduced byte for byte.

The deposited records are those of the first run, which the paper uses.

## Times

On a recent Mac with 8 processes:
- `alg_diagrams.py` and `validate_plumbing.py`: about 30 s each;
- `crosscheck_u.py` and `cw_general.py`: about 1 min each.

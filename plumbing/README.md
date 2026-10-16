# Double branched covers, Neumann normal forms and correction terms (Sections 2.2 and 3)

Run the scripts from this folder, and those in `neumann/` from `neumann/`. Scripts marked (Sage) run with
`sage -python` (see the top-level README). Each script writes the JSON file of the same name next to it.

## 1. The double branched cover and its recognition

$\Sigma(K)$ is built as the 2-fold cyclic cover of the orbifold obtained by $(2,0)$-filling the exterior of $K$ in
SnapPy. Its triangulation is then recognised with Regina.

| Script | Output | Content |
| --- | --- | --- |
| `reach_orb.py` | `reach_u12_orb.json` | the 437 knots with previous range $[1,2]$: Regina name of $\Sigma(K)$ from 4 randomised triangulations. `reach_u12.json` is the earlier filled-cover route of `reach.py`, kept for comparison |
| `reach_u23_orb.py` | `reach_u23_orb.json` | the 666 knots with previous range $[2,3]$ or $[2,4]$, other than the 47 non-alternating knots with $u=u_a=3$ (those are handled by `classify47.py`) |
| `recognise_more.py` | `recognise_more.json` | 60 more randomisations for the covers not recognised yet |
| `fresh_targets.py` | `recognise_more_fresh.json` | the covers newly recognised as graph manifolds; input of `neumann/run_fresh.py` |
| `recognise_retry.py` | `recognise_retry.json` | up to 300 randomisations with `simplifyExhaustive`, also for the mirror |
| `verify_hyp.py` (Sage) | `verify_hyp.json` | verified hyperbolicity. A short geodesic is drilled and refilled, which gives another triangulation of $\Sigma(K)$ with a positively oriented solution; that solution is then certified with interval arithmetic |
| `jsj.py` | `jsj.json` | incompressible tori among the vertex normal surfaces (quadrilateral coordinates), and the pieces cut along them |
| `jsj_verify.py` (Sage) | `jsj_verify.json` | verified hyperbolicity of the pieces, which gives positive simplicial volume |
| `pachner_recognise.py` | `pachner_recognise.json` | recognition by a search of the Pachner graph (height at most 3) |

## 2. Normal forms and obstructions (`neumann/`)

| Script | Output | Content |
| --- | --- | --- |
| `formula.py` | | the Regina name of the boundary of a plumbing tree, computed by formula |
| `calibrate.py` | `calibrate.json` | checks `formula.py` against Regina's recognition of triangulated plumbings (random pair and triple trees) |
| `invert.py` | | building blocks of normal forms: Hirzebruch–Jung legs, chains, node weights from $\lvert\det\rvert=D$ |
| `precheck.py` | `precheck_u12.json`, `precheck_u23.json` | the normal forms of $\Sigma(K)$ and $-\Sigma(K)$ from the Regina names; verdict `APPLICABLE` (one is negative definite) or `NOT APPLICABLE` (both have $p^+\ge1$) |
| `precheck_sigma0.py` | `precheck_sigma0.json` | the same for the 32 graph-manifold covers with $\sigma=0$ that no other run processed |
| `run_invert.py` | | module: orientation from the linking pairing, and the $u=1$ test (Ni–Wu formula; graded roots and the mapping cone for one bad vertex) |
| `run_fresh.py` | `run_fresh.json` | normal form and tests for `recognise_more_fresh.json`. For $[1,2]$ the $u=1$ test; for $[2,3]$ and $[2,4]$ Owens's test (`../owens_u3_graph.py`) |
| `run_retry13.py`, `run_pachner4.py` | `run_retry13.json`, `run_pachner4.json` | the same for the covers recognised by `recognise_retry.py` and `pachner_recognise.py` |
| `rederive_all.py` | `rederive_all.json` | the 43 new values of Sections 3.2 and 3.3 from the normal form alone: verdict, Casson–Walker identity, and an independent identification by triangulating the plumbing |
| `normal_case3.py` | `normal_case3.json` | the same for the 96 knots of class (3) whose plumbing was first found by the searches of Section 3 below |
| `verify_fresh.py` | `verify_fresh.json` | checks of `run_fresh.json`: independent Ni–Wu test, independent Owens step for non-cyclic groups (13n2956), Traczyk's criterion from $V_K(e^{i\pi/3})$ in exact arithmetic |

## 3. Earlier plumbing searches and verification

| Script | Output | Content |
| --- | --- | --- |
| `star_u1.py`, `star_u1_orb.py` | `star_u1.json`, `star_u1_orb.json` | star-shaped plumbings found by search for Seifert fibred covers, and the $u=1$ test |
| `graph_u1.py` | `graph_u1.json` | the same for graph manifolds with two or three Seifert pieces; `NO PLUMBING` when the search finds none |
| `owens_u3_graph.py` | `owens_u3_targets.json`, `owens_u3_controls.json` | Owens's test with plumbing correction terms. `controls` runs the 469 Montesinos knots of Section 3.4 |
| `owens_u3_verify.py` | `owens_u3_verify.json` | independent Owens check, cyclic groups |
| `owens_general_group.py` | `owens_general_group.json` | independent Owens check for non-cyclic groups: every form, exact class minima, every group isomorphism (13n2933) |
| `crosscheck.py` | `crosscheck.json`, `crosscheck_controls.json` | recognition by the orbifold route; Casson–Walker identity with $\sum d+2\sum\operatorname{rank}HF_{\mathrm{red}}$ |
| `graph_verify.py` | `graph_verify.json` | calibration of the closed Casson–Walker formula $L(\Gamma)$ on star plumbings and Montesinos controls |

The verdicts used in the paper come from the normal forms of Section 2, as `../paper/make_tables.py` asserts. The
records of this section supply the Regina names used as input.

## 4. The 47 non-alternating knots with $u=u_a=3$

`classify47.py` applies the same functions to the 47 knots, one stage at a time:
- `recognise`
- `hyp` (Sage)
- `jsj`
- `piece` (Sage)
- `normal` (normal forms; for the two knots with a negative-definite normal form, Owens's test)
- `verify` (triangulation identification and the independent Owens check)
- `summary`

| Class | Knots |
| --- | --- |
| hyperbolic | 39 |
| toroidal with a hyperbolic piece | 5 |
| no negative-definite plumbing tree | 1 (13n1662) |
| negative-definite normal form | 2 (13n3004, 13n3005) |

For 13n3004 and 13n3005, Owens's test excludes $u=2$ in both orientations.

## Modules

- `general.py`: names, legs and linking forms.
- `plumbing.py`: plumbing matrices, the plumbing link, triangulated plumbed manifolds.
- `d_owens.py`: correction terms of a negative-definite plumbing tree, $d=\max(c^2+\lvert V\rvert)/4$.

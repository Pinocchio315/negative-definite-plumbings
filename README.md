# Negative-definite plumbings and new unknotting numbers

This repository contains the code and the computational records for the paper
**Negative-Definite Plumbings and New Unknotting Numbers** by Seong-Jin Lee.

The paper continues [arXiv:2609.09861][v2] (*Computation of unknotting numbers: which knot breaks the
Bernhard–Jablan conjecture?*). Its code is in [unknotting-tait-graphs][v2code]. Unknotting ranges are those of the
final table of arXiv:2609.09861v2, and knot data come from KnotInfo through `database_knotinfo` 2026.8.1.

## Results

**Double branched covers (Section 3.1).**
- The study covers the 1,150 non-alternating knots with at most thirteen crossings whose unknotting number lay in
  $[1,2]$, $[2,3]$ or $[2,4]$.
- The second column below also counts the 47 non-alternating knots with $u=u_a=3$ (Section 4.1).

| $\Sigma(K)$ | $u\in[1,2]$ | $u\in[2,3],[2,4]$ | total |
| --- | ---: | ---: | ---: |
| hyperbolic (verified) | 187 | 303 | 490 |
| toroidal, with a hyperbolic JSJ piece | 69 | 80 | 149 |
| graph manifold; neither $\Sigma(K)$ nor $-\Sigma(K)$ bounds a negative-definite plumbing tree | 102 | 197 | 299 |
| $\Sigma(K)$ or $-\Sigma(K)$ bounds a negative-definite plumbing tree | 59 | 121 | 180 |
| toroidal, JSJ pieces not identified | 20 | 12 | 32 |
| **total** | **437** | **713** | **1,150** |

For 938 of these knots, no orientation of $\Sigma(K)$ bounds a negative-definite plumbing tree. A graph manifold falls
in the third or fourth row according to its Neumann normal forms, by Neumann's theorem on the positive index of
plumbing trees.

**New values.**

| Result | Knots | Method |
| --- | --- | --- |
| $u=2$ (§3.2) | 41 knots, previous range $[1,2]$ | normal-form correction terms exclude $u=1$: Ni–Wu formula (35 knots), reduced Floer homology and the mapping cone (6) |
| $u=3$ (§3.3) | 13n2933, 13n2956 | Owens's $u=2$ obstruction with normal-form correction terms, and Traczyk's criterion |
| $u=u_a=3$ (§4.1) | 66 knots (19 alternating, 47 non-alternating), previous range $[2,3]$ | rank-two linking pairing obstruction (Borodzik–Friedl) and three crossing changes |
| $u=5$ (§4.2) | 13a4770, 13a4787, 13a4801 | theorem of McCoy–Park and a finite check of one diagram |
| second explicit counterexample to the Bernhard–Jablan conjecture (§5) | 13n1587 | $u=2$, and every crossing change in a minimal diagram gives a knot with $u\ge2$ |

In total, 112 unknotting numbers are new. For 13n3004 and 13n3005, two of the 66 knots, Owens's test also excludes
$u=2$. This gives a second proof for them.

## Layout

```text
.
├── data/                 inputs: ranges of arXiv:2609.09861v2, the new values, KnotInfo extracts
├── lib/                  shared modules (from the arXiv:2609.09861 repository, byte-identical) and helpers
├── plumbing/             Sections 2.2 and 3: Sigma_2, hyperbolicity, JSJ, Neumann normal forms, obstructions
│   └── neumann/          normal forms of Sigma_2 and -Sigma_2, tests, re-derivation and verification
├── rank_two/             Section 4.1: the rank-two linking pairing test, 66 knots, controls
│   └── alternating_audit/  independent implementation for the 19 alternating knots
├── special_alternating/  Section 4.2: u = 5 for three special alternating knots
├── bernhard_jablan/      Section 5: 13n1587 (minimal diagrams, neighbours, upper bound)
└── paper/                make_tables.py and the generated LaTeX tables
```

Each folder has a README that lists its scripts and records. Most JSON files are deposited outputs of the scripts
next to them. Running a script rewrites its output in place, so work in a copy when comparing.

## Environment

The computations used Python 3.13.7 and SageMath 10.7 on macOS.

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

- **Python packages.** The pinned versions are in `requirements.txt`: SnapPy 3.3.2, Spherogram 2.4.1, Regina 7.4.1,
  knot_floer_homology 1.2.2, NumPy 2.5.2, SymPy 1.14.0, NetworkX 3.6.1 and `database_knotinfo` 2026.8.1.
- **KnotInfo data.** The KnotInfo CSV of `database_knotinfo` 2026.8.1 has SHA-256
  `e37cfb46060b0ba5165a730e482275c8a82460c6a462653e662188cfe952598d`.
- **SageMath.** Verified hyperbolicity uses SnapPy's interval arithmetic inside SageMath. This concerns
  `plumbing/verify_hyp.py`, `plumbing/jsj_verify.py`, and the `hyp` and `piece` stages of `plumbing/classify47.py`.
  - Run them with `sage -python`.
  - If SnapPy, Regina and `database_knotinfo` are installed only in the virtual environment, first set
    `VENV_SITE_PACKAGES=$PWD/.venv/lib/python3.13/site-packages`.
- **Signatures.** Signatures follow KnotInfo, where the right-handed trefoil has $\sigma=-2$.

## Regenerating the tables of the paper

```sh
python -B paper/make_tables.py
```

The script reads the deposited records and writes `paper/tables/*.tex`. One of these is `numbers.tex`, which holds
macros for every count in the text. It takes a few seconds and asserts the conditions the paper relies on:
- Every new value comes from an obstruction computed with the normal form, satisfies the Casson–Walker identity, and
  was identified by an independent triangulation.
- Every knot of class (3) of Section 3.1 has a negative-definite normal form, with at most two bad vertices.
- For every knot of class (2), two different trees are accepted, the normal forms of $\Sigma(K)$ and $-\Sigma(K)$,
  and neither is negative definite.
- The classification is complete: no knot is left unclassified.
- The rank-two test has no false positive among the negative controls and detects all positive controls.

## Re-running the computations

Run every script from its own folder, as in the folder READMEs. The table lists the times of the re-runs on a
recent Mac with 8–11 processes.

| Step | Command | Time |
| --- | --- | --- |
| classification of the double branched covers of the 47 non-alternating knots with $u=u_a=3$ | `cd plumbing; python -B classify47.py recognise`, `sage -python classify47.py hyp`, `python -B classify47.py jsj`, `sage -python classify47.py piece`, `python -B classify47.py normal`, `python -B classify47.py verify`, `python -B classify47.py summary` | about 5 min |
| normal forms deciding class (2) | `cd plumbing/neumann; python -B precheck.py u12; python -B precheck.py u23` | 2 s per knot |
| both normal forms for each knot of class (2) | `cd plumbing/neumann; python -B both_orientations.py` | about 4 min |
| the 43 new values from the normal form alone | `cd plumbing/neumann; python -B rederive_all.py` | about 1 min |
| the 96 knots of class (3) first found by plumbing searches | `cd plumbing/neumann; python -B normal_case3.py` | about 2 min |
| rank-two test on the 66 knots, controls | `cd rank_two; python -B check66.py; python -B controls_bf2.py` | 20 s |
| $u=5$ | `cd special_alternating; python -B verify3.py` | 75 s |
| neighbours of 13n1587 | `cd bernhard_jablan/code; python -B neighbor_lower.py ../results/bjset_13n_1587.json` | seconds |
| upper bound for 13n1587 | `cd bernhard_jablan/upper_bound; python -B check.py` | 3.5 min |

Every row above was re-run in this layout, and each reproduced the deposited records. The only differences were
the following:
- run times;
- tuples written as lists;
- for 13n1662, the randomised recognition named the cover at a later attempt, with the same name.

## Paper and code

| Paper | Claim | Code | Records |
| --- | --- | --- | --- |
| §2.2 | the plumbing graph ↔ Regina notation; normal forms from Regina names | `plumbing/neumann/formula.py`, `calibrate.py`, `precheck.py` | `calibrate.json` |
| §3.1 | recognition of $\Sigma(K)$ | `plumbing/reach_orb.py`, `reach_u23_orb.py`, `recognise_more.py`, `recognise_retry.py` | `reach_*_orb.json`, `recognise_*.json` |
| §3.1 | verified hyperbolicity | `plumbing/verify_hyp.py` (Sage) | `verify_hyp.json` |
| §3.1 | incompressible tori and hyperbolic pieces | `plumbing/jsj.py`, `jsj_verify.py` (Sage), `pachner_recognise.py` | `jsj.json`, `jsj_verify.json`, `pachner_recognise.json` |
| §3.1 | normal forms of $\pm\Sigma(K)$ | `plumbing/neumann/precheck.py`, `precheck_sigma0.py`, `run_fresh.py`, `run_retry13.py`, `run_pachner4.py`, `normal_case3.py`, `both_orientations.py` | `plumbing/neumann/*.json` |
| §3.1 | the 47 non-alternating knots of §4.1 | `plumbing/classify47.py` | `classify47_*.json` |
| §3.2 | $u=1$ excluded: Ni–Wu, mapping cone | `plumbing/neumann/run_invert.py` (`u1_test`), `lib/lower_bounds/montesinos/` | `rederive_all.json` |
| §3.3 | $u=2$ excluded: Owens, Traczyk | `plumbing/owens_u3_graph.py` (`run_one`), `lib/lower_bounds/owens/` | `rederive_all.json`, `classify47_normal.json` |
| §3.4 | identification, Casson–Walker identity | `plumbing/neumann/rederive_all.py`, `plumbing/owens_u3_graph.py` (`L_closed`), `plumbing/graph_verify.py` | `rederive_all.json`, `graph_verify.json` |
| §3.4 | independent implementations | `bernhard_jablan/code/u1_independent.py`, `owens_independent.py`, `plumbing/owens_general_group.py`, `plumbing/neumann/verify_fresh.py` | `owens_general_group.json`, `verify_fresh.json`, `classify47_verify.json` |
| §3.4 | controls of Owens's test | `plumbing/owens_u3_graph.py controls` | `owens_u3_controls.json` |
| §4.1 | rank-two linking pairing test | `rank_two/bf2.py`, `extend.py`, `robust47.py`, `upper47.py`, `check66.py`, `controls_bf2.py` | `rank_two/*.json` |
| §4.1 | second implementation for the alternating knots | `rank_two/alternating_audit/independent_audit/code/independent_hnf.py` | `hnf_results.json` |
| §4.2 | $u=5$ | `special_alternating/verify3.py` | `verify3.json` |
| §5 | minimal diagrams of 13n1587 and 12n563 | `bernhard_jablan/code/enum_diagrams.py`, `neighbors.py`, `verify_completeness.py`, `cover_invariant.py` | `bernhard_jablan/results/` |
| §5 | neighbours | `neighbor_lower.py`, `composite_neighbors.py`, `cert_identify*.py`, `neighbor_rigor.py`, `u1_independent.py` | `bernhard_jablan/results/` |
| §5 | $u(13n1587)\le2$ | `bernhard_jablan/upper_bound/check.py`; the braid certificate of Brittenham–Hermiller: `upper_bounds/verify_braid_certificate.py` of the arXiv:2609.09861 repository | `check.json`, `hits.jsonl` |
| tables | all tables and counts | `paper/make_tables.py` | `paper/tables/` |

## Provenance

- **The single method of the paper.** Plumbings were first found by searches over candidate trees, recorded in
  `plumbing/star_u1*.json`, `graph_u1.json` and `owens_u3_targets.json`.
  - All verdicts of the paper come from the Neumann normal form:
    - the 43 new values: `rederive_all.json`;
    - the knots of class (3): `run_fresh.json`, `run_retry13.json`, `run_pachner4.json`, `normal_case3.json`,
      `precheck_sigma0.json`;
    - the knots of class (2): `precheck_*.json`, with both normal forms listed in `both_orientations.json`;
    - the 47 non-alternating knots with $u=u_a=3$: `classify47_*.json`.
  - `make_tables.py` asserts this.
  - The search records still supply Regina names, and they are the inputs of the normal-form runs.
- **Regenerated records.** `plumbing/neumann/precheck_u12.json` and `precheck_u23.json` were regenerated with the
  final `precheck.py`, which also tries mixed orientations of the Seifert pieces. All 189 verdicts and minimal
  positive indices are the same as in the first run.
- **Rebuilt input.** `plumbing/recognise_more_fresh.json` was first written by a one-off command.
  `plumbing/fresh_targets.py` rebuilds it with identical content.
- **Copied modules.**
  - `lib/lower_bounds/` and `lib/tait_graphs/tait/` are byte-identical copies of the modules of the
    arXiv:2609.09861 repository.
  - `data/greene_sweep_2026-09-08.jsonl.gz` is its Greene-model sweep.
  - `data/arxiv_2609.09861v2_main.tex` is the source of arXiv:2609.09861v2, which `rank_two/parse_paper.py` reads.
- **Starting ranges.** `lib/ranges.py` contains the function `load_ranges` of the project module that fixed the
  starting ranges, copied verbatim.
- **Rank-two scripts.** The scripts in `rank_two/` were first run on 5 October 2026. Their records are from an
  unchanged re-run on 15 October 2026.
- **Upper bound for 13n1587.** The 15-crossing diagram of 13n1587 in `bernhard_jablan/upper_bound/hits.jsonl` was
  found by a randomised Reidemeister search, which is not included. `check.py` verifies these certificates
  independently of the search.

## Citation

Seong-Jin Lee, *Negative-Definite Plumbings and New Unknotting Numbers* (2026).

The preceding paper: Seong-Jin Lee, *Computation of unknotting numbers: which knot breaks the Bernhard–Jablan
conjecture?*, [arXiv:2609.09861][v2] [math.GT] (2026).

[v2]: https://arxiv.org/abs/2609.09861
[v2code]: https://github.com/Pinocchio315/unknotting-tait-graphs

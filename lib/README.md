# Shared modules

`lower_bounds/` and `tait_graphs/tait/` are byte-identical copies of the corresponding modules of the
[arXiv:2609.09861 repository](https://github.com/Pinocchio315/unknotting-tait-graphs). There they are documented and
validated.

| Module | Used for |
| --- | --- |
| `lower_bounds/owens/owens_obstruction.py` | class minima $m_Q$, the forms $\widetilde Q$ of Owens's theorem, matching of correction terms |
| `lower_bounds/owens/traczyk.py`, `kinfo.py` | Traczyk's criterion from $V_K(e^{i\pi/3})$; KnotInfo access |
| `lower_bounds/montesinos/montesinos_u1.py` | correction terms of star plumbings, the Ni–Wu test of $u=1$ over all affine identifications |
| `lower_bounds/montesinos/hf_red.py` | ranks of $HF_{\mathrm{red}}$ from the graded root, for at most one bad vertex |
| `lower_bounds/linking_pairing.py` | linking pairings of finite abelian groups |
| `lower_bounds/owens/taitgraph.py`, `_compat.py`, `tait_graphs/tait/` | Goeritz matrices of alternating diagrams (controls of Owens's test) |

Helpers written in this project:

| Module | Content |
| --- | --- |
| `ranges.py` | `load_ranges`: the starting ranges (the final table of arXiv:2609.09861v2 and the values of this project) |
| `exact_helpers.py` | crossing changes on PD codes and exact determinants (from `../bernhard_jablan/code/verify_completeness.py`) |
| `owens_traczyk_8_6_4_1.py` | used only for `qt_ind` and `gram_reduce`: the matrices $\widetilde Q$ and their Gram reduction for the independent Owens checks |

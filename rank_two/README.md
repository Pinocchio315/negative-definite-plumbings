# Rank-two linking pairings: knots with $u=u_a=3$ (Section 4.1)

Run the scripts from this folder.

**The test.** `bf2.py` decides whether a knot can be turned into a knot with Alexander polynomial one by two crossing
changes. This is the obstruction of Borodzik and Friedl, in the form of Lemma 4.1 of the paper.
- $C$ runs over the binary forms of BF Lemma 5.2, with $\lvert\det C\rvert=\det K$ and both signs of $\det C$.
- Signature: $\operatorname{sign}C=\sigma(K)+e_1+e_2$.
- Diagonal residues: modulo 4 they are $e_1,e_2$, with an even off-diagonal entry, tested over the 96 elements of
  $GL(2,\mathbb Z/4)$.
- Linking pairing: $2\,l(K)(v_i,v_j)\equiv(C^{-1})_{ij}$ for some $v_1,v_2\in H_1(\Sigma(K))$.
- The test uses the signed mode with reduced forms: `analyse(V, mode='signed', reduced=True)`.

| Script | Output | Content |
| --- | --- | --- |
| `extend.py` | `extend.jsonl` | the test on the 2,264 knots whose KnotInfo range is $[1,3]$, $[1,4]$, $[2,3]$, $[2,4]$ or $[2,5]$; 96 are obstructed |
| `parse_paper.py` | `paper_results.json` | the exact values and improved ranges of arXiv:2609.09861v2, read from `../data/arxiv_2609.09861v2_main.tex` |
| `select47.py` | `new47.txt` | the obstructed 13-crossing non-alternating knots that are not exact in v2 and have upper bound 3 |
| `robust47.py` | `robust47.json` | the test again, with the KnotInfo Seifert matrix and with one computed from the PD code; also the test without sign conditions |
| `upper47.py` | `upper47.json` | three crossing changes in the KnotInfo diagram that give the unknot (knot Floer homology of total rank one, $\pi_1\cong\mathbb Z$); none with two or fewer |
| `check66.py` | `check66.json` | `robust47.work` and `upper47.work` on all 66 knots of `../data/new66_u3.txt` |
| `controls_bf2.py` | `controls_bf2.json` | negative controls: 9,673 knots with $u\le2$ (v2) or $u_a\le2$ (KnotInfo), none excluded. Positive controls: the 25 knots of BF §6.2 and their mirrors, all excluded |

`published_negative_controls.json` lists the 25 knots of BF §6.2 with their KnotInfo Seifert matrices.

**Results.** For all 66 knots, the test excludes two crossing changes, so $u_a\ge3$:
- 57 are excluded without the sign conditions.
- The other 9 have $\lvert\sigma\rvert=4$; for them the signature condition fixes $e_1=e_2$, and the residue and
  linking conditions then exclude every form.

Three crossing changes give the unknot in each case. The 96 obstructed knots of `extend.jsonl` are the 66 knots and
30 knots that already had $u\ge3$; none of them contradicts an upper bound.

**Second implementation (`alternating_audit/`).** `independent_audit/code/independent_hnf.py` treats the 19
alternating knots and their mirrors without importing any of the code above. It uses:
- integer Hermite normal forms for $\operatorname{coker}(V+V^T)$;
- exact rational arithmetic for the pairing;
- its own enumeration of forms and of $GL(2,\mathbb Z/4)$.

Its inputs are the frozen cohort files in `inputs/`, `bf/` and `derived/`.

```sh
cd alternating_audit
BF_AUDIT_OUTPUT_DIR=$(mktemp -d) python -B independent_audit/code/independent_hnf.py
```

The re-run reproduces every verdict of `independent_audit/results/hnf_summary.json`. Only the recorded code hash
changes.

The scripts were first run on 5 October 2026, and the deposited outputs are from an unchanged re-run on
15 October 2026. `check66.py` and `controls_bf2.py` take about 10 seconds each.

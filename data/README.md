# Inputs

| File | Content |
| --- | --- |
| `v2_final_u_table.json` | the final unknotting ranges of arXiv:2609.09861v2 for the prime knots with at most 13 crossings (KnotInfo 2026.8.1 with the 104 DKT-adopted upper bounds rolled back, plus the results of v2) |
| `arxiv_2609.09861v2_results.json` | the exact values and improved ranges of arXiv:2609.09861v2 (equal to `../rank_two/paper_results.json`) |
| `arxiv_2609.09861v2_main.tex` | the LaTeX source of arXiv:2609.09861v2, read by `../rank_two/parse_paper.py` |
| `knotinfo_rollback_list_2026-09-07.csv` | the 104 rolled-back KnotInfo entries |
| `knotinfo_2026-09-16_unknotting.json` | KnotInfo unknotting data of 16 September 2026, used for novelty checks |
| `knotinfo_2026-10-05_new_values.json` | the `unknotting_number` entries of KnotInfo version 2026.10.5 for the 112 knots of Theorem 1.2; all are ranges. Written by `knotinfo_2026_10_05.py` from `database_knotinfo/csv_data/knotinfo_data_complete.csv` of the wheel `database_knotinfo-2026.10.5-py3-none-any.whl` (PyPI, downloaded on 15 October 2026), which is not included. The file records the SHA-256 of the wheel and of the CSV; the script checks the latter against the RECORD file of the wheel |
| `new66_u3.txt` | the 66 knots with $u=u_a=3$ (Section 4.1) |
| `new_values_20261015.json` | the 43 values of Sections 3.2 and 3.3 (`strict=true`) and the three values $u=5$ (`strict=false`) |
| `greene_sweep_2026-09-08.jsonl.gz` | Greene-model correction terms from the arXiv:2609.09861 repository (`results/greene/sweep_2026-09-08.jsonl.gz`) |

Knot diagrams, Seifert matrices, signatures, determinants and Jones polynomials come from the package
`database_knotinfo==2026.8.1`.

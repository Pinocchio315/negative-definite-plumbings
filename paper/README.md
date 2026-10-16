# Tables of the paper

`make_tables.py` reads the deposited records and writes the LaTeX fragments in `tables/`:

| File | Content |
| --- | --- |
| `sigma2_map.tex` | the classification of the double branched covers (Section 3.1, Table 1) |
| `new_u2.tex`, `new_u3.tex` | the new values $u=2$ and $u=3$ (Sections 3.2 and 3.3), with $\Sigma(K)$ and the size of the normal form |
| `u3_66.tex` | the 66 knots with $u=u_a=3$ (Section 4.1) |
| `nbr_13n1587.tex` | the knots obtained by one crossing change in a minimal diagram of 13n1587 |
| `numbers.tex` | LaTeX macros for the counts used in the text |

Run from the repository root:

```sh
python -B paper/make_tables.py
```

No count is typed by hand: every number in the text comes from `numbers.tex`.

# Special alternating knots with $u=5$ (Section 4.2)

`verify3.py` writes `verify3.json`. It treats 13a4770, 13a4787 and 13a4801, which have $\sigma=-8$ and previous range
$[4,5]$:
- It checks that all thirteen crossings of the reduced alternating KnotInfo diagram have the same sign.
- It checks that no set of at most four crossing changes gives the unknot: the determinant is computed exactly, then
  the total rank of knot Floer homology is checked.
- It records a set of five crossing changes that gives the unknot: knot Floer homology of total rank one and
  $\pi_1\cong\mathbb Z$.

With the theorem of McCoy and Park (arXiv:2603.10732, Thm 1), this gives $u=5$. Runtime: about 75 seconds.

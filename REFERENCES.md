# References

This curated bibliography supports the public v0.1 software: exact candidate
generation, specialization scoring, Sage/PARI/eclib analysis, and proof-status
language.

## High-Rank Elliptic-Curve Search

### Mestre (1982)

Jean-François Mestre.  
**Construction d'une courbe elliptique de rang >= 12.**  
*C. R. Acad. Sci. Paris Ser. I Math.*, 295 (1982), 643-644.

Used here for historical background on high-rank elliptic-curve construction
and finite-prime search heuristics.

### Nagao (1992)

Koh-ichi Nagao.  
**Examples of elliptic curves over Q with rank >= 17.**  
*Proceedings of the Japan Academy, Series A, Mathematical Sciences*, 68(9)
(1992), 287-289.  
DOI: 10.3792/pjaa.68.287.

Used here for the high-rank specialization-search lineage of the
Mestre-Nagao scoring method.

### Nagao (1994)

Koh-ichi Nagao.  
**Construction of high-rank elliptic curves.**  
*Kobe Journal of Mathematics*, 11(2) (1994), 211-219.

Used here for high-rank elliptic-curve construction and computational-search
background.

## Specialization Background

### Silverman (1983)

Joseph H. Silverman.  
**Heights and the specialization map for families of abelian varieties.**  
*Journal für die reine und angewandte Mathematik*, 342 (1983), 197-211.

Used here for background on specialization of Mordell-Weil groups in families.

## Published Benchmark Family

### Elkies and Klagsbrun (2020)

Noam D. Elkies and Zev Klagsbrun.  
**New Rank Records For Elliptic Curves Having Rational Torsion.**  
*Open Book Series*, 4 (2020), 233-250.  
arXiv:2003.00077.  
DOI: 10.2140/obs.2020.4.233.  
https://arxiv.org/abs/2003.00077

Used here for:

- the explicit two-parameter elliptic K3 family with rational 2-torsion used
  in the public specialization benchmark;
- the specialization at `u = 2/5`;
- the published rank-19 benchmark at `t = 11860/97527`;
- the nine generic Mordell-Weil sections used in the benchmark workflow;
- the Mestre-Nagao score
  `S(t,B) = sum_{p < B, good reduction} log(N_p(E_t)/p)`;
- local precomputation as a function of `t mod p`;
- fixed-denominator sieving and staged filtering.

The public implementation is independent project code built from the published
mathematics and documented algorithms; it is not copied from the authors'
source code.

### Elkies 2026 — rank-17 elliptic K3 fibration

**Source:** Noam D. Elkies, *An elliptic K3 surface $X/\mathbb{Q}(t)$ with Mordell-Weil rank 17, I: Formulas for $X$ and base changes of ranks 18 and 19*, arXiv:2608.25406v1, 26 August 2026.

**Type:** Primary source.

**Claims used in this repository:**

- Gives an explicit elliptic K3 surface over $\mathbb{Q}(t)$ with 17 independent sections.
- Provides the explicit Weierstrass model and the 17 section coordinates used by the rank-17 reproduction.
- Gives the height-pairing calculation and the published $17\times17$ Gram matrix with determinant $948$.
- Identifies
  $$
  t=-9529/5471
  $$
  as the specialization yielding Elkies's published rank-at-least-28 curve.
- Describes explicit base changes producing generic ranks 18 and 19.

**Repository connection:**  
The scripts under `examples/elkies_2026_rank17/` reproduce the published rank-17 fibration and the published rank-at-least-28 specialization. This repository uses the paper as a source for reproduction of published results and reproducibility infrastructure, not for unpublished specialization searches or new rank claims.

## Computational Tools

### SageMath

The Sage Developers.  
**SageMath, the Sage Mathematics Software System.**  
https://www.sagemath.org/

Used here for elliptic curves, rational arithmetic, canonical heights, point
construction, and interfaces to arithmetic libraries. A version-specific
citation should be added for archived public releases when the release
environment is finalized.

### PARI/GP

The PARI Group.  
**PARI/GP.**  
Université de Bordeaux.  
https://pari.math.u-bordeaux.fr/

Used here through SageMath for arithmetic diagnostics, including rank-related
computations where available. PARI diagnostics are not treated as rank
certificates without the corresponding mathematical justification.

### eclib / mwrank

John Cremona.  
**eclib / mwrank.**  
https://github.com/JohnCremona/eclib

Used here through SageMath for Mordell-Weil and descent computations. Release
documentation should cite the exact eclib version used when reporting
reproducible rank-bound computations.

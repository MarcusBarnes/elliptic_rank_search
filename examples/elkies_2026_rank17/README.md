# Reproducing Elkies's 2026 rank-17 elliptic K3 fibration

This example reproduces a **published result** and is intended for the public repository. It deliberately stops at reproduction and verification; it contains no unpublished high-rank specialization searches, candidate lists, record-oriented heuristics, or new lattice-theoretic claims.

## Source

Noam D. Elkies, *An elliptic K3 surface X/Q(t) with Mordell-Weil rank 17, I: Formulas for X and base changes of ranks 18 and 19*, arXiv:2608.25406v1, 26 August 2026.

Relevant material:

- Theorem 4 / equation (1): the Weierstrass model `y^2 = x^3 - 27 S(t) x + (27/4) T(t)`.
- Section 2.1: the 17 polynomial `x`-coordinates.
- Equation (2): signs of the leading coefficients of the corresponding `y`-coordinates.
- Lemmas 5 and 6: the height calculation for integral sections.
- Theorem 4: the published 17-by-17 height Gram matrix, with determinant 948.

## What the script verifies

`elkies_2026_rank17.sage`:

1. constructs the published polynomials `S(t)` and `T(t)`;
2. constructs all 17 published `x_i(t)`;
3. obtains each `y_i(t)` as an exact polynomial square root of the curve equation, choosing the sign published by Elkies;
4. checks the fully printed `(x_1,y_1)` pair against the paper;
5. verifies all 17 sections satisfy the Weierstrass equation exactly;
6. recomputes every height pairing from Elkies's Lemma 6;
7. verifies the resulting Gram matrix equals the published matrix;
8. verifies that its determinant is `948`, proving the 17 sections are independent.

Run from the repository root with Sage available:

```bash
sage examples/elkies_2026_rank17/elkies_2026_rank17.sage
```

Expected final output includes:

```text
17 published sections verified exactly
Gram matrix matches published matrix
determinant = 948
sections are independent
```

## Research boundary

This reproduction is public infrastructure. Extensions intended to discover new rank records, unusually small conductors/heights, unpublished base changes, improved specialization heuristics, or new mathematical structure belong in the private research workflow until independently verified and deliberately released.


## Published rank-28 specialization

The companion script

```bash
sage examples/elkies_2026_rank17/elkies_2026_rank28_specialization.sage
```

reproduces the published rank-at-least-28 fiber.

Elkies (2026, Section 2.3) identifies

\[
t=-9529/5471
\]

as the fiber yielding his 2006 rank-at-least-28 curve. The script specializes
the rank-17 fibration at this value and checks **exactly over \(\mathbf Q\)**
that the resulting elliptic curve is isomorphic to the published model

\[
y^2+xy+y=x^3-x^2+a_4x+a_6,
\]

with the coefficients printed by Klagsbrun--Sherman--Weigandt (2016).

It also:

- specializes the 17 generic sections and checks they remain rational points;
- checks the 28 published rational points on the published rank-28 model;
- recomputes their Néron--Tate height-pairing matrix and numerical regulator
  as a reproducibility check of the published independence claim.

The exact-rank statement is **not unconditional**: Klagsbrun, Sherman, and
Weigandt prove that the curve has rank exactly 28 subject to GRH for number
fields. The unconditional published statement needed here is rank at least 28,
coming from the 28 independent rational points.

This benchmark is intentionally reproduction-only. Searching the fibration
for new high-rank fibers or other record features belongs to the private
research workflow.

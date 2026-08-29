# Small Reproducible Workflows

These workflows are intended for the future public repository. They avoid large
searches and keep exact arithmetic separate from Sage/PARI diagnostics.

## Workflow 1: Pure-Python Integral Smoke Search

Command:

```bash
python3 search_integral.py \
  --umin -1 --umax 1 \
  --vmin -3 --vmax 3 \
  --min-points 3 \
  --top 3 \
  --output small_integral_demo.json
```

Expected behavior:

- writes `data/candidates/small_integral_demo.json`;
- reports at least one nonsingular candidate curve;
- the top candidate in the current implementation has `A=-115`, `B=417`,
  discriminant `18776000`, `distinct_x_count=3`, and `visible_point_count=5`.

Regression-testable exact invariants:

- every listed visible point satisfies `y^2 + x*y = x^3 + A*x + B`;
- the generalized discriminant is nonzero;
- the points arise from the exact `q(U,V)` collinearity relation.

## Workflow 2: Pure-Python Rational Auxiliary Point Check

Command:

```bash
python3 -m unittest discover -s tests
```

Expected behavior:

- all tests pass without SageMath;
- tests check that `Z=1` rational auxiliary points agree with the integral
  formula;
- tests check exact rational point and discriminant identities.

Regression-testable exact invariants:

- `search_rational.make_aux_point(U,V,1).q == search_integral.q_value(U,V)`;
- `q = y^2 + x*y - x^3` for rational auxiliary points;
- duplicate `Z` values are canonicalized by `parse_z_values`.

## Workflow 3: Published Specialization Benchmark

This workflow requires SageMath.

Command:

```bash
sage search_specializations.py \
  --u 2/5 \
  --t 11860/97527 \
  --output data/benchmarks/benchmark_u2_5_t11860_97527.json
```

Expected behavior:

- writes a small JSON benchmark artifact;
- constructs the published Elkies-Klagsbrun specialization at `u=2/5`,
  `t=11860/97527`;
- recovers the published generic-section x-coordinates where exact square-root
  recovery succeeds.

Regression-testable exact invariants:

- `5-u^2` is a rational square;
- the constructed curve is nonsingular;
- recovered section coordinates satisfy the specialized curve equation exactly.

Do not treat the recovered generic sections or numerical height diagnostics as a
record-rank certificate. Rank certification requires separate rigorous evidence.

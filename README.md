# Elliptic Rank Search

[![DOI](https://zenodo.org/badge/1350985562.svg)](https://doi.org/10.5281/zenodo.22167424)

Research software for reproducible computational searches for elliptic curves over
$\mathbb{Q}$ with large Mordell-Weil rank.

The primary public objective of this project is to develop computational methods
for finding promising elliptic curves over $\mathbb{Q}$, analyzing them with
exact and heuristic tools, and clearly separating candidate generation from
rigorous rank certification.

This repository does not claim a new rank record.

## Mathematical Scope

For an elliptic curve $E/\mathbb{Q}$, the Mordell-Weil group $E(\mathbb{Q})$
is finitely generated. Its rank is one of the central arithmetic invariants of
the curve, and finding curves with unusually large rank is a long-running
computational problem.

This software supports several distinct stages of that process:

- constructing candidate curves with visible rational points;
- verifying those points using exact arithmetic;
- using local scores and Sage/PARI diagnostics to prioritize candidates;
- analyzing candidate Mordell-Weil subgroups;
- preparing data that may later support rigorous rank certification.

These stages should not be conflated. A list of visible rational points may be
linearly dependent. A PARI rank estimate, numerical height-pairing rank,
heuristic Mestre-Nagao score, or Selmer upper bound is not by itself a certified
Mordell-Weil rank lower bound.

## What the Software Currently Does

The initial public release focuses on stable, reproducible Track-A tooling:

- exact pure-Python search for curves
  `y^2 + x*y = x^3 + A*x + B` with visible integral points from a structured
  residue-class ansatz;
- exact pure-Python extension to rational points with denominator parameter `Z`;
- SageMath analysis of generated candidate curves;
- specialization construction and NumPy Mestre-Nagao scoring for a published
  Elkies-Klagsbrun family benchmark;
- small reproducibility tests and examples.

Additional unpublished exploratory components are intentionally excluded from
the initial public release unless and until they are mature enough for
publication.

## Reproduced From Literature

The specialization workflow is based on published methods of Mestre, Nagao, and
Elkies-Klagsbrun.

In particular, the repository contains a small benchmark for the
Elkies-Klagsbrun specialization at

```text
u = 2/5
t = 11860/97527
```

This benchmark is included as a reproducibility and provenance artifact. It is
not presented as a new mathematical result.

Any public statement about current record ranks should be treated as a
time-sensitive claim and verified against citable sources before release.

## Quick Start

Run a small pure-Python exact search:

```bash
python3 search_integral.py \
  --umin -1 --umax 1 \
  --vmin -3 --vmax 3 \
  --min-points 3 \
  --top 3 \
  --output small_integral_demo.json
```

Expected current top candidate:

```text
A=-115
B=417
discriminant=18776000
distinct_x_count=3
visible_point_count=5
```

Run the lightweight test suite:

```bash
python3 -m unittest discover -s tests
```

Analyze an included exact example with SageMath:

```bash
sage sage_analyze.py example_curve.json
```

Reproduce the published specialization benchmark with SageMath:

```bash
sage search_specializations.py \
  --u 2/5 \
  --t 11860/97527 \
  --output data/benchmarks/benchmark_u2_5_t11860_97527.json
```

See `docs/reproducibility/SMALL_WORKFLOWS.md` for additional reproducibility
notes and expected outputs.

## Environment

The pure-Python baseline requires Python 3.11 or newer.

Specialization scoring requires NumPy.

Arithmetic analysis requires SageMath with PARI/GP and eclib/mwrank available.
Depending on platform, SageMath may be installed through conda-forge, a system
package manager, or the official SageMath distribution. See `environment.yml`
for a starting conda environment.

## Architecture

The public workflow is organized around a separation between search, analysis,
and certification:

```text
candidate construction
        |
        v
exact rational-point verification
        |
        v
cheap arithmetic / local scoring
        |
        v
candidate shortlist
        |
        v
SageMath arithmetic analysis
        |
        +--> visible rational points
        +--> height-pairing diagnostics
        +--> PARI diagnostics
        +--> eclib / mwrank bounds
        |
        v
candidates requiring rigorous certification
```

The early stages are intentionally small and exact where practical. Expensive
or heuristic searches should produce artifacts recording parameters, seeds,
source references, and proof-status metadata.

## Reproducibility

Small workflows are documented in
`docs/reproducibility/SMALL_WORKFLOWS.md`.

The default tests avoid SageMath so that they can run quickly in ordinary
Python environments. Sage-dependent workflows are documented separately because
standard CI environments may not provide a practical Sage/eclib installation.

For computational claims, this project aims to preserve enough information to
distinguish:

- exact computations;
- reproducible heuristic searches;
- numerical or diagnostic evidence;
- conditional mathematical statements;
- rigorously certified results.

## Limitations

- The baseline collinearity search is exact but not intended to scale to the
  largest searches of interest.
- Search results are candidates, not rank certificates.
- The Mestre-Nagao score is a heuristic ranking signal.
- Sage/PARI/mwrank outputs must be interpreted according to the method used:
  lower bounds, upper bounds, exact results, and conditional claims are
  mathematically different.
- Large generated search campaigns and unpublished candidate sets are
  intentionally excluded from the initial public release.

## References

Core references include work by Mestre, Nagao, Silverman, and
Elkies-Klagsbrun. See `REFERENCES.md` for the curated public bibliography.

## AI and Computational Tool Disclosure

This project uses artificial intelligence as part of an interactive
computational research workflow.

OpenAI's ChatGPT has been used as a research and software-engineering assistant
for activities including literature exploration, mathematical discussion,
hypothesis generation, algorithm and experiment design, code generation and
debugging, interpretation of computational results, and preparation of
technical documentation.

The research questions, project direction, selection and execution of
experiments, evaluation of mathematical significance, and decisions about which
results and claims to retain are the responsibility of the human researcher.
Computational outputs and AI-generated suggestions are not treated as
mathematical proof merely because they were produced by an automated system.
Claims intended to carry mathematical weight are expected to be supported
independently through exact computation, reproducible evidence, formal
verification, proof, or appropriate reference to the literature.

This disclosure is made in the spirit of the **Leiden Declaration on Artificial
Intelligence and Mathematics**, particularly its recommendations concerning
transparent disclosure of automated tools, human responsibility for correctness,
appropriate attribution, and the distinction between automated assistance and
human authorship.

AI systems are not listed as authors and are not assigned responsibility for
the mathematical results in this repository.

## Citation

Use `CITATION.cff` for citation metadata. Tagged public releases are intended
to be archived through Zenodo so that citable, versioned software records can be
preserved. Users should cite the exact software version they used.

## License

This public release is distributed under the BSD-3-Clause license. See
`LICENSE`.

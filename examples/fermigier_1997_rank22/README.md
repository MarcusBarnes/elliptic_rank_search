# Fermigier 1997 rank ≥ 22 reproduction

This example reproduces selected computational features of Stéfane Fermigier's published 1997 rank ≥ 22 elliptic curve construction.

The purpose of the example is **reproducibility and provenance**. It is not a new rank claim and it is not the project's unpublished specialization-search code.

## Source

Stéfane Fermigier, *Une courbe elliptique définie sur $\mathbb{Q}$ de rang $\ge 22$*, Acta Arithmetica **82**(4) (1997), 359–363. DOI: `10.4064/aa-82-4-359-363`.

## What is reproduced

For the selected split sextic with roots

$$
0,55,314,378,1007,1036,
$$

the script uses Fermigier's published specialization parameter

$$
t_{\mathrm{paper}}=\frac{19754}{39}.
$$

The symmetric-shift implementation in this repository uses

$$
T=2t_{\mathrm{paper}}=\frac{39508}{39}
$$

in

$$
q_6(x-T)q_6(x+T)=g(x)^2-r(x).
$$

The script checks exactly that $\deg r=4$, verifies the twelve forced rational $x$-coordinates $a_i\pm T$ on the quartic $y^2=r(x)$, constructs the Jacobian, and verifies that it is $\mathbb{Q}$-isomorphic to Fermigier's published curve

$$
y^2+xy+y=x^3-940299517776391362903023121165864x
+10707363070719743033425295515449274534651125011362.
$$

It also reproduces Fermigier's historical Mestre–Nagao score checkpoints using the convention that $M$ denotes the first $M$ primes, with the $p=2$ contribution omitted:

| $M$ | published score |
|---:|---:|
| 50 | 29.49 |
| 100 | 44.12 |
| 200 | 57.54 |
| 400 | 81.51 |
| 1000 | 105.17 |
| 2000 | 122.76 |

## Run

```bash
sage examples/fermigier_1997_rank22/reproduce_rank22.sage
```

A successful run writes `reproduction_result.json` beside the script and exits with status 0.

## Proof status

Fermigier's paper proves rank at least $22$. This example does **not** independently re-certify all 22 independent points. It reproduces the construction, the published specialization up to $\mathbb{Q}$-isomorphism, and the historical score calibration. Those distinctions are intentional.

No unpublished candidate parameters or large specialization-search results are included here.

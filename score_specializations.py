#!/usr/bin/env python3
"""
Staged Mestre--Nagao specialization scorer for the Elkies--Klagsbrun
Z/2Z-torsion elliptic K3 family.

MATHEMATICAL PROVENANCE
=======================

The search strategy implemented here is based primarily on:

[EK2020]
    N. D. Elkies and Z. Klagsbrun,
    "New Rank Records For Elliptic Curves Having Rational Torsion",
    Open Book Series 4 (2020), 233--250.
    arXiv:2003.00077
    DOI: 10.2140/obs.2020.4.233

In particular:

* EK2020, Section 2 defines the score

      S(t,B) = sum_{p < B, good reduction} log(N_p(E_t)/p).

* EK2020, Section 3 observes (following Nagao) that a_p(E_t) depends
  only on t mod p, so local values can be precomputed.

* EK2020, Section 3 describes staged screening with increasing prime
  bounds B_0 <= ... <= B_m and increasingly selective cutoffs.

* EK2020, Section 3.1 describes the fixed-denominator sieve used here:
  for t=a/b, fix b and an interval of numerators a, and add periodic
  local score arrays into one counter array.

* EK2020, Section 3.1 stores rounded values
      round(D * log(N_p(E_t)/p))
  with D=1024 so the sieve can use integer addition efficiently.
  We follow that fixed-point idea, but use NumPy int32 counters by
  default for safety and portability rather than assuming 16-bit
  bounds.

Historical sources cited by EK2020 include:

[Mestre1982]
    J.-F. Mestre,
    "Construction de courbes elliptiques sur Q de rang >= 12",
    C. R. Acad. Sci. Paris Ser. I Math. 295 (1982), 643--644.

[Nagao1992]
    K.-I. Nagao,
    "Examples of elliptic curves over Q with rank >= 17",
    Proc. Japan Acad. Ser. A Math. Sci. 68 (1992), 287--289.

[Nagao1994]
    K.-I. Nagao,
    "Construction of high-rank elliptic curves",
    Kobe J. Math. 11 (1994), 211--219.

This implementation is an independent, open-source reproduction/adaptation
of the published computational ideas. It is not source code from EK2020.

IMPLEMENTATION STATUS
=====================

This version optimizes the *search-region scoring*:

1. Precompute local score tables s_p(t mod p).
2. Sieve an entire fixed-denominator numerator interval at once using
   periodic NumPy array additions.
3. Keep only the strongest candidates after the first stage.
4. Extend the score of survivors at larger prime bounds by table lookup.

The local table generator still counts points for each t in F_p using
straightforward finite-field arithmetic. That remains the main scaling
bottleneck at large B and is intentionally isolated so that it can later
be replaced by a faster native/vectorized implementation without changing
the search logic.

Example benchmark
=================

Known rank-19 specialization from EK2020:

    u = 2/5
    t = 11860/97527

Example calibration:

    python3 score_specializations_optimized.py \
        --u 2/5 \
        --denominator 97527 \
        --amin 8000 \
        --amax 16000 \
        --stage-bounds 300,1000 \
        --stage-keep 1000,100 \
        --benchmark-numerator 11860 \
        --output specialization_scores_staged.json
"""

from __future__ import annotations

import argparse
import json
import math
import time
from fractions import Fraction
from math import gcd
from pathlib import Path

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit(
        "This optimized scorer requires NumPy. "
        "Install it in the active environment before running."
    ) from exc


REFERENCE_METADATA = [
    {
        "key": "EK2020",
        "authors": "Noam D. Elkies; Zev Klagsbrun",
        "title": "New Rank Records For Elliptic Curves Having Rational Torsion",
        "year": 2020,
        "journal": "Open Book Series",
        "volume": "4",
        "pages": "233-250",
        "arxiv": "2003.00077",
        "doi": "10.2140/obs.2020.4.233",
        "used_for": [
            "Mestre-Nagao score",
            "local t mod p precomputation",
            "staged prime-bound filtering",
            "fixed-denominator sieving",
            "D=1024 fixed-point score representation",
        ],
    },
    {
        "key": "Mestre1982",
        "authors": "Jean-Francois Mestre",
        "title": "Construction de courbes elliptiques sur Q de rang >= 12",
        "year": 1982,
        "journal": "C. R. Acad. Sci. Paris Ser. I Math.",
        "volume": "295",
        "pages": "643-644",
        "used_for": ["historical origin of high-rank finite-prime scoring"],
    },
    {
        "key": "Nagao1992",
        "authors": "Koh-Ichi Nagao",
        "title": "Examples of elliptic curves over Q with rank >= 17",
        "year": 1992,
        "journal": "Proc. Japan Acad. Ser. A Math. Sci.",
        "volume": "68",
        "pages": "287-289",
        "used_for": ["family specialization version of the Mestre-Nagao method"],
    },
    {
        "key": "Nagao1994",
        "authors": "Koh-Ichi Nagao",
        "title": "Construction of high-rank elliptic curves",
        "year": 1994,
        "journal": "Kobe J. Math.",
        "volume": "11",
        "pages": "211-219",
        "used_for": ["high-rank elliptic-curve construction/search background"],
    },
]


def parse_fraction(text: str) -> Fraction:
    return Fraction(text)


def frac_str(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def resolve_output_path(filename: str, default_dir: str) -> Path:
    path = Path(filename)
    if path.parent == Path("."):
        path = Path(default_dir) / path
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def parse_int_list(text: str):
    vals = [int(x.strip()) for x in text.split(",") if x.strip()]
    if not vals:
        raise argparse.ArgumentTypeError("Expected a comma-separated list of integers.")
    return vals


def primes_below(n: int):
    if n <= 2:
        return []

    sieve = bytearray(b"\x01") * n
    sieve[0:2] = b"\x00\x00"

    for p in range(2, int(n**0.5) + 1):
        if sieve[p]:
            start = p * p
            sieve[start:n:p] = b"\x00" * (((n - 1 - start) // p) + 1)

    return [p for p in range(2, n) if sieve[p]]


def inv_mod(a: int, p: int) -> int:
    return pow(a % p, -1, p)


def fraction_mod_p(x: Fraction, p: int):
    den = x.denominator % p
    if den == 0:
        return None
    return (x.numerator % p) * inv_mod(den, p) % p


def A_mod_p(u: int, t: int, p: int) -> int:
    u2 = u*u % p
    u3 = u2*u % p
    u4 = u2*u2 % p
    u5 = u4*u % p
    u6 = u3*u3 % p
    u7 = u6*u % p
    u8 = u4*u4 % p

    c4 = (u8 - 18*u6 + 163*u4 - 1152*u2 + 4096) % p
    c3 = (3*u7 - 35*u5 - 120*u3 + 1536*u) % p
    c2 = (u8 - 13*u6 + 32*u4 - 152*u2 + 1536) % p
    c1 = (u7 + 3*u5 - 156*u3 + 672*u) % p
    c0 = (3*u6 - 33*u4 + 112*u2 - 80) % p

    t2 = t*t % p
    t3 = t2*t % p
    t4 = t2*t2 % p

    return (c4*t4 + c3*t3 + c2*t2 + c1*t + c0) % p


def B_factors_mod_p(u: int, t: int, p: int):
    def b1(uu, tt):
        return ((uu*uu + uu - 8)*tt + (-uu + 2)) % p

    def b3(uu, tt):
        return ((uu*uu - uu - 8)*tt + (uu*uu + uu - 10)) % p

    def b5(uu, tt):
        return ((uu*uu - 7*uu + 8)*tt + (-uu*uu + uu + 2)) % p

    def b7(uu, tt):
        return ((uu*uu + 5*uu + 8)*tt + (uu*uu + 3*uu + 2)) % p

    nu = (-u) % p
    nt = (-t) % p

    return [
        b1(u, t),
        (-b1(nu, nt)) % p,
        b3(u, t),
        (-b3(nu, nt)) % p,
        b5(u, t),
        (-b5(nu, nt)) % p,
        b7(u, t),
        (-b7(nu, nt)) % p,
    ]


def B_mod_p(u: int, t: int, p: int) -> int:
    value = 1
    for factor in B_factors_mod_p(u, t, p):
        value = value * factor % p
    return value


def discriminant_mod_p(a2: int, a4: int, p: int) -> int:
    return (16 * a4 * a4 * (a2*a2 - 4*a4)) % p


def legendre_table(p: int):
    table = np.full(p, -1, dtype=np.int8)
    table[0] = 0
    y = np.arange(1, p, dtype=np.int64)
    table[(y*y) % p] = 1
    return table


def point_count(p: int, A: int, B: int, chi) -> int:
    """
    Count E(F_p) for
        y^2 = x^3 + 2*A*x^2 + B*x.

    This is the deliberately simple local-table backend. Search-region
    sieving is optimized; local table construction remains the next major
    optimization target.
    """
    x = np.arange(p, dtype=np.int64)
    rhs = (x*x % p * x + (2*A % p)*(x*x % p) + B*x) % p
    return int(p + 1 + chi[rhs].sum(dtype=np.int64))


def build_local_table(
    u: Fraction,
    p: int,
    denominator: int,
    fixed_scale: int,
):
    """
    Build one local table indexed by t mod p.

    Values are fixed-point integers
        round(D * log(N_p(E_t)/p)).
    Bad-reduction entries contribute zero, matching omission from S(t,B).
    """
    u_mod = fraction_mod_p(u, p)

    if u_mod is None:
        return None, "u_denominator"

    if denominator % p == 0:
        return None, "search_denominator"

    chi = legendre_table(p)
    scores = np.zeros(p, dtype=np.int32)
    good = np.zeros(p, dtype=np.bool_)

    for t_mod in range(p):
        A = A_mod_p(u_mod, t_mod, p)
        B = B_mod_p(u_mod, t_mod, p)
        a2 = 2*A % p

        if discriminant_mod_p(a2, B, p) == 0:
            continue

        Np = point_count(p, A, B, chi)
        scores[t_mod] = int(round(fixed_scale * math.log(Np / p)))
        good[t_mod] = True

    return {
        "p": p,
        "b_inv": inv_mod(denominator, p),
        "scores": scores,
        "good": good,
    }, None


def periodic_update(table, amin: int):
    """
    Construct the p-periodic fixed-point update pattern for numerators
    a = amin+i.

    This is the mathematical periodicity exploited in EK2020 §3.1.
    """
    p = table["p"]
    b_inv = table["b_inv"]

    residues = np.arange(p, dtype=np.int64)
    t_residues = ((amin % p) + residues) * b_inv % p
    return table["scores"][t_residues]


def add_periodic(scores, period):
    """
    Add one periodic local-score pattern into the full counter array.

    NumPy performs the bulk additions in native code. This is a portable
    approximation to the vectorized counter-array sieve described in EK2020.
    """
    n = len(scores)
    p = len(period)

    full = n // p
    rem = n % p

    if full:
        scores[:full*p].reshape(full, p)[:] += period

    if rem:
        scores[full*p:] += period[:rem]


def reduced_mask(amin: int, amax: int, denominator: int):
    nums = np.arange(amin, amax + 1, dtype=np.int64)
    return np.fromiter(
        (gcd(int(a), denominator) == 1 for a in nums),
        dtype=np.bool_,
        count=len(nums),
    )


def score_one_from_tables(a: int, tables):
    total = 0
    good_count = 0
    bad_count = 0

    for table in tables:
        p = table["p"]
        t_mod = (a % p) * table["b_inv"] % p

        if table["good"][t_mod]:
            total += int(table["scores"][t_mod])
            good_count += 1
        else:
            bad_count += 1

    return total, good_count, bad_count


def parse_args():
    ap = argparse.ArgumentParser()

    ap.add_argument("--u", default="2/5")
    ap.add_argument("--denominator", type=int, default=97527)
    ap.add_argument("--amin", type=int, default=8000)
    ap.add_argument("--amax", type=int, default=16000)

    ap.add_argument(
        "--stage-bounds",
        type=parse_int_list,
        default=parse_int_list("300,1000"),
        help=(
            "Increasing prime bounds, e.g. 1000,2048,8192. "
            "The first stage uses the fixed-denominator sieve; later stages "
            "extend only surviving candidates."
        ),
    )

    ap.add_argument(
        "--stage-keep",
        type=parse_int_list,
        default=parse_int_list("1000,100"),
        help=(
            "Candidates retained after each stage. Must have the same "
            "number of entries as --stage-bounds."
        ),
    )

    ap.add_argument(
        "--fixed-scale",
        type=int,
        default=1024,
        help="Fixed-point denominator D. EK2020 used D=1024.",
    )

    ap.add_argument(
        "--benchmark-numerator",
        type=int,
        default=11860,
    )

    ap.add_argument(
        "--output",
        default="specialization_scores_staged.json",
    )

    return ap.parse_args()


def main():
    args = ap = parse_args()

    bounds = args.stage_bounds
    keeps = args.stage_keep

    if len(bounds) != len(keeps):
        raise SystemExit("--stage-bounds and --stage-keep must have equal length")

    if any(b <= 3 for b in bounds):
        raise SystemExit("All stage bounds must be > 3")

    if any(bounds[i] >= bounds[i+1] for i in range(len(bounds)-1)):
        raise SystemExit("--stage-bounds must be strictly increasing")

    if any(k <= 0 for k in keeps):
        raise SystemExit("--stage-keep values must be positive")

    if any(keeps[i] < keeps[i+1] for i in range(len(keeps)-1)):
        raise SystemExit("--stage-keep must be nonincreasing")

    if args.amin > args.amax:
        raise SystemExit("--amin must be <= --amax")

    if args.denominator <= 0:
        raise SystemExit("--denominator must be positive")

    u = parse_fraction(args.u)
    b = args.denominator
    D = args.fixed_scale

    all_primes = [p for p in primes_below(bounds[-1]) if p != 2]

    print("Staged Elkies--Klagsbrun Mestre--Nagao scorer")
    print(f"u = {frac_str(u)}")
    print(f"fixed denominator b = {b}")
    print(f"numerator interval = [{args.amin}, {args.amax}]")
    print(f"stage bounds = {bounds}")
    print(f"stage keep counts = {keeps}")
    print(f"fixed-point scale D = {D}")
    print()

    mask = reduced_mask(args.amin, args.amax, b)
    numerators = np.arange(args.amin, args.amax + 1, dtype=np.int64)
    reduced_indices = np.flatnonzero(mask)

    print(f"reduced fractions in interval = {len(reduced_indices)}")
    print()

    tables = []
    skipped = {
        "u_denominator": [],
        "search_denominator": [],
    }

    next_prime_index = 0
    previous_bound = 2

    # Candidate state after first stage is represented by absolute array index.
    survivor_indices = reduced_indices.copy()
    survivor_fixed_scores = None

    stage_results = []

    for stage_no, (bound, keep) in enumerate(zip(bounds, keeps), start=1):
        started = time.time()

        print(f"Stage {stage_no}: extending local tables to p < {bound}", flush=True)

        new_tables = []

        while next_prime_index < len(all_primes):
            p = all_primes[next_prime_index]

            if p >= bound:
                break

            next_prime_index += 1

            if p < previous_bound:
                continue

            table, reason = build_local_table(
                u=u,
                p=p,
                denominator=b,
                fixed_scale=D,
            )

            if table is None:
                skipped[reason].append(p)
                continue

            tables.append(table)
            new_tables.append(table)

            if len(new_tables) % 25 == 0:
                print(
                    f"  built {len(new_tables)} new local tables "
                    f"(latest p={p})",
                    flush=True,
                )

        if stage_no == 1:
            # EK2020 §3.1-style fixed-denominator sieve.
            full_scores = np.zeros(
                args.amax - args.amin + 1,
                dtype=np.int32,
            )

            for table in new_tables:
                period = periodic_update(table, args.amin)
                add_periodic(full_scores, period)

            reduced_scores = full_scores[survivor_indices]

            order = np.argsort(reduced_scores)[::-1]
            retain = min(keep, len(order))

            survivor_indices = survivor_indices[order[:retain]]
            survivor_fixed_scores = reduced_scores[order[:retain]].astype(
                np.int64
            )

        else:
            # EK2020 staged trick: only extend survivors at later bounds.
            increments = np.zeros(len(survivor_indices), dtype=np.int64)

            for j, idx in enumerate(survivor_indices):
                a = int(numerators[idx])
                subtotal = 0

                for table in new_tables:
                    p = table["p"]
                    t_mod = (a % p) * table["b_inv"] % p

                    if table["good"][t_mod]:
                        subtotal += int(table["scores"][t_mod])

                increments[j] = subtotal

            survivor_fixed_scores += increments

            order = np.argsort(survivor_fixed_scores)[::-1]
            retain = min(keep, len(order))

            survivor_indices = survivor_indices[order[:retain]]
            survivor_fixed_scores = survivor_fixed_scores[order[:retain]]

        elapsed = time.time() - started

        stage_entry = {
            "stage": stage_no,
            "prime_bound": bound,
            "new_local_table_count": len(new_tables),
            "cumulative_local_table_count": len(tables),
            "survivors": len(survivor_indices),
            "elapsed_seconds": elapsed,
        }

        stage_results.append(stage_entry)

        print(
            f"  survivors = {len(survivor_indices)}; "
            f"new tables = {len(new_tables)}; "
            f"time = {elapsed:.2f}s"
        )
        print()

        previous_bound = bound

    rows = []

    for idx, fixed_score in zip(survivor_indices, survivor_fixed_scores):
        a = int(numerators[idx])

        # Recompute good/bad counts only for final survivors.
        _, good, bad = score_one_from_tables(a, tables)

        rows.append({
            "numerator": a,
            "denominator": b,
            "t": frac_str(Fraction(a, b)),
            "fixed_score": int(fixed_score),
            "score": float(fixed_score) / D,
            "good_prime_count": good,
            "bad_prime_count": bad,
        })

    rows.sort(key=lambda r: r["fixed_score"], reverse=True)

    benchmark = None

    if args.benchmark_numerator:
        a = args.benchmark_numerator

        if gcd(a, b) == 1:
            fixed_score, good, bad = score_one_from_tables(a, tables)

            # Rank benchmark within the final candidate list when present;
            # otherwise report how many final survivors exceed it.
            better_final = sum(
                1 for row in rows
                if row["fixed_score"] > fixed_score
            )

            benchmark = {
                "known_rank": 19,
                "numerator": a,
                "denominator": b,
                "t": frac_str(Fraction(a, b)),
                "fixed_score": fixed_score,
                "score": fixed_score / D,
                "good_prime_count": good,
                "bad_prime_count": bad,
                "final_survivors_with_higher_score": better_final,
                "survived_final_stage": any(
                    row["numerator"] == a for row in rows
                ),
            }

    output_path = resolve_output_path(
        args.output,
        "data/candidates",
    )

    payload = {
        "provenance": {
            "references": REFERENCE_METADATA,
            "notes": [
                "Search score and sieve/staging design adapted from EK2020 Sections 2, 3, and 6.",
                "This is an independent open-source implementation.",
                "Fixed-point scale D=1024 follows the value reported in EK2020 Section 3.1.",
            ],
        },
        "method": {
            "score": "sum_{p<B, good} log(N_p(E_t)/p)",
            "fixed_point_scale": D,
            "first_stage": "fixed-denominator periodic sieve",
            "later_stages": "direct local-table lookup on survivors",
            "local_table_backend": "NumPy finite-field point counting",
        },
        "search": {
            "u": frac_str(u),
            "denominator": b,
            "amin": args.amin,
            "amax": args.amax,
            "stage_bounds": bounds,
            "stage_keep": keeps,
            "initial_reduced_fraction_count": int(len(reduced_indices)),
            "skipped_primes": skipped,
        },
        "stages": stage_results,
        "benchmark": benchmark,
        "candidates": rows,
    }

    with output_path.open("w") as f:
        json.dump(payload, f, indent=2)

    if benchmark is not None:
        print("Known rank-19 benchmark at final bound:")
        print(f"  t = {benchmark['t']}")
        print(f"  score ~= {benchmark['score']:.12f}")
        print(f"  good primes = {benchmark['good_prime_count']}")
        print(f"  bad primes = {benchmark['bad_prime_count']}")
        print(
            f"  survived final stage = "
            f"{benchmark['survived_final_stage']}"
        )
        print(
            "  final survivors with higher score = "
            f"{benchmark['final_survivors_with_higher_score']}"
        )
        print()

    print("Top final candidates:")

    for i, row in enumerate(rows[:20], start=1):
        marker = ""
        if (
            args.benchmark_numerator
            and row["numerator"] == args.benchmark_numerator
        ):
            marker = "  <-- known rank-19 benchmark"

        print(
            f"{i:>3}. t={row['t']} "
            f"score~={row['score']:.12f} "
            f"good={row['good_prime_count']} "
            f"bad={row['bad_prime_count']}"
            f"{marker}"
        )

    print()
    print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()

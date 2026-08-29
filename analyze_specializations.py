#!/usr/bin/env sage
"""
Analyze high-scoring Elkies--Klagsbrun specializations.

This is the arithmetic-analysis stage following the Mestre--Nagao scorer.

For each selected specialization, the script can:

1. Construct the exact specialized curve.
2. Recover the 9 published generic sections.
3. Compute their canonical-height span/regulator.
4. Optionally run Sage point_search().
5. Optionally run PARI rank diagnostics.
6. Optionally run eclib/mwrank in Selmer-only mode on the global minimal
   model, obtaining a rigorous rank upper bound without global Selmer-class searches.
7. Compare that upper bound with a generic-rank reference and a campaign
   target rank, without treating numerical section span as a proof.
8. Optionally run Sage saturation on recovered points only when explicitly requested.

METHODOLOGICAL PROVENANCE
=========================

Primary family/search reference:

    N. D. Elkies and Z. Klagsbrun,
    "New Rank Records For Elliptic Curves Having Rational Torsion",
    Open Book Series 4 (2020), 233--250.
    arXiv:2003.00077
    DOI: 10.2140/obs.2020.4.233

The direct Mordell--Weil/descent interface used here is SageMath's interface
to John Cremona's eclib library (mwrank).

Important interpretation:
    * numerical height span is a screening diagnostic, not a proof by itself;
    * PARI rank output is recorded as a diagnostic and should not be used alone
      as a certified Mordell-Weil rank statement;
    * mwrank rank() is a lower bound in general;
    * mwrank rank_bound() is an upper bound obtained from 2-descent;
    * mwrank certain() reports whether the rank has been determined.

This script is independent project tooling built around published mathematics
and open-source arithmetic software.
"""

from __future__ import annotations

import argparse
import json
import time
from fractions import Fraction
from math import isqrt
from pathlib import Path

from sage.all import QQ, RR, EllipticCurve
from sage.libs.eclib.interface import mwrank_EllipticCurve


# ---------------------------------------------------------------------------
# Basic utilities
# ---------------------------------------------------------------------------

def Q(s):
    if isinstance(s, Fraction):
        return s
    return Fraction(str(s))


def frac_str(x: Fraction) -> str:
    return (
        str(x.numerator)
        if x.denominator == 1
        else f"{x.numerator}/{x.denominator}"
    )


def resolve_output_path(filename: str, default_dir: str) -> Path:
    path = Path(filename)

    if path.parent == Path("."):
        path = Path(default_dir) / path

    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def rational_square_root(x: Fraction):
    if x < 0:
        return None

    rn = isqrt(x.numerator)
    rd = isqrt(x.denominator)

    if rn * rn == x.numerator and rd * rd == x.denominator:
        return Fraction(rn, rd)

    return None


# ---------------------------------------------------------------------------
# Elkies--Klagsbrun family
# ---------------------------------------------------------------------------

def A_polynomial(u: Fraction, t: Fraction) -> Fraction:
    c4 = u**8 - 18*u**6 + 163*u**4 - 1152*u**2 + 4096
    c3 = 3*u**7 - 35*u**5 - 120*u**3 + 1536*u
    c2 = u**8 - 13*u**6 + 32*u**4 - 152*u**2 + 1536
    c1 = u**7 + 3*u**5 - 156*u**3 + 672*u
    c0 = 3*u**6 - 33*u**4 + 112*u**2 - 80

    return c4*t**4 + c3*t**3 + c2*t**2 + c1*t + c0


def B1(u, t):
    return (u**2 + u - 8)*t + (-u + 2)


def B3(u, t):
    return (u**2 - u - 8)*t + (u**2 + u - 10)


def B5(u, t):
    return (u**2 - 7*u + 8)*t + (-u**2 + u + 2)


def B7(u, t):
    return (u**2 + 5*u + 8)*t + (u**2 + 3*u + 2)


def B_factors(u, t):
    return [
        B1(u, t),
        -B1(-u, -t),
        B3(u, t),
        -B3(-u, -t),
        B5(u, t),
        -B5(-u, -t),
        B7(u, t),
        -B7(-u, -t),
    ]


def B_product(u, t):
    value = Fraction(1)

    for factor in B_factors(u, t):
        value *= factor

    return value


def compatible_m_values(u, sqrt_s):
    a = u - 2
    b = Fraction(2)
    c = u + 2

    roots = []

    if a == 0:
        if b != 0:
            roots = [-c / b]
    else:
        disc = b*b - 4*a*c
        sd = rational_square_root(disc)

        if sd is not None:
            roots = [
                (-b + sd)/(2*a),
                (-b - sd)/(2*a),
            ]

    out = []

    for m in roots:
        den = m*m + 1

        if den == 0:
            continue

        u_check = 2*(m*m - m - 1)/den
        s_check = (m*m + 4*m - 1)/den

        if u_check == u and abs(s_check) == abs(sqrt_s):
            out.append((m, s_check))

    return out


def published_section_xs(u, t, m):
    b1,b2,b3,b4,b5,b6,b7,b8 = B_factors(u, t)

    return [
        -b1*b2*b3*b6,
        -b1*b2*b4*b5,
        4*b1*b2*b5*b6,
        b1*b3*b4*b6,
        -b1*b3*b4*b7,
        b1*b3*b4*b8,
        b1*b3*b5*b6,
        -b1*b5*b6*b7,
        -(m - 1)**2 * b1*b2*b3*b8,
    ]


# ---------------------------------------------------------------------------
# Generic sections and numerical span
# ---------------------------------------------------------------------------

def recover_generic_sections(E, u, t, A, B):
    square_term = Fraction(5) - u*u
    sqrt_s = rational_square_root(square_term)

    if sqrt_s is None:
        return {
            "success": False,
            "error": "5-u^2 is not a rational square",
            "selected_m": None,
            "points": [],
            "point_rows": [],
            "failed_sections": [],
        }

    m_pairs = compatible_m_values(u, sqrt_s)

    if not m_pairs:
        return {
            "success": False,
            "error": "No compatible rational m found",
            "selected_m": None,
            "points": [],
            "point_rows": [],
            "failed_sections": [],
        }

    m_pairs.sort(key=lambda pair: pair[1] != sqrt_s)
    m = m_pairs[0][0]

    xs = published_section_xs(u, t, m)

    points = []
    rows = []
    failures = []

    for i, xf in enumerate(xs, start=1):
        rhs = xf**3 + 2*A*xf**2 + B*xf
        yf = rational_square_root(rhs)

        if yf is None:
            failures.append({
                "section": i,
                "x": frac_str(xf),
                "error": "RHS is not a rational square",
            })
            continue

        xq = QQ(xf.numerator) / QQ(xf.denominator)
        yq = QQ(yf.numerator) / QQ(yf.denominator)

        try:
            P = E(xq, yq)
        except Exception as exc:
            failures.append({
                "section": i,
                "x": frac_str(xf),
                "error": repr(exc),
            })
            continue

        points.append(P)
        rows.append({
            "section": i,
            "x": str(P[0]),
            "y": str(P[1]),
        })

    return {
        "success": len(points) == 9,
        "error": None if len(points) == 9 else f"Recovered {len(points)}/9 sections",
        "selected_m": frac_str(m),
        "points": points,
        "point_rows": rows,
        "failed_sections": failures,
    }


def numerical_span(E, points, rel_tol=1e-12):
    if not points:
        return {
            "estimated_span": 0,
            "height_determinant": None,
            "regulator_of_points": None,
            "height_eigenvalues": [],
        }

    H = E.height_pairing_matrix(points)
    eigs = sorted([abs(RR(x)) for x in H.eigenvalues()], reverse=True)

    scale = max(eigs) if eigs else RR(0)
    tol = scale * RR(rel_tol) if scale else RR(rel_tol)

    span = sum(1 for x in eigs if x > tol)

    regulator = None

    try:
        regulator = E.regulator_of_points(points)
    except Exception:
        pass

    return {
        "estimated_span": int(span),
        "height_determinant": str(H.det()),
        "regulator_of_points": None if regulator is None else str(regulator),
        "height_eigenvalues": [str(x) for x in eigs],
    }


def unique_points(points):
    chosen = []
    keys = set()

    for P in points:
        if P.is_zero():
            continue

        Qneg = -P

        key1 = (str(P[0]), str(P[1]))
        key2 = (str(Qneg[0]), str(Qneg[1]))
        canonical = min(key1, key2)

        if canonical in keys:
            continue

        keys.add(canonical)
        chosen.append(P)

    return chosen


def extra_points_only(generic_points, found_points):
    generic_keys = set()

    for P in generic_points:
        generic_keys.add((P[0], P[1]))
        Qneg = -P
        generic_keys.add((Qneg[0], Qneg[1]))

    extras = []

    for P in found_points:
        if P.is_zero():
            continue

        if (P[0], P[1]) in generic_keys:
            continue

        extras.append(P)

    return unique_points(extras)


# ---------------------------------------------------------------------------
# Direct mwrank/eclib analysis
# ---------------------------------------------------------------------------

def run_mwrank_screen(
    Emin,
    generic_rank_reference,
    target_rank,
    first_limit,
    second_limit,
    second_descent,
    verbose,
):
    """
    Fast Selmer-only upper-bound screen.

    IMPORTANT:
        The numerical span of specialized published sections is NOT used as a
        rigorous lower bound.  This routine only records the rigorous mwrank
        upper bound and compares it with:
          * the published generic-rank reference for the family; and
          * the comparison target rank.

    Thus:
        rank_jump_upper = upper_bound - generic_rank_reference

    is only an upper bound on the possible specialization jump.  It is useful
    for triage and population-level analysis, but does not prove that the jump
    occurs.
    """
    started = time.time()

    result = {
        "enabled": True,
        "mode": "selmer_only_upper_bound_screen",
        "generic_rank_reference": int(generic_rank_reference),
        "target_rank": int(target_rank),
        "first_limit": first_limit,
        "second_limit": second_limit,
        "second_descent": second_descent,
        "rank_upper_bound": None,
        "rank_jump_upper": None,
        "rank_jump_candidate": None,
        "can_reach_target": None,
        "target_rejected_by_upper_bound": None,
        "exact_rank_certified": False,
        "minimal_model_ainvs": [str(x) for x in Emin.a_invariants()],
        "error": None,
        "elapsed_seconds": None,
        "note": (
            "Upper-bound screening only. Numerical section-span estimates are "
            "stored separately and are not used as rank proofs."
        ),
    }

    try:
        ainvs = [int(x) for x in Emin.a_invariants()]
        M = mwrank_EllipticCurve(ainvs, verbose=verbose)

        M.two_descent(
            verbose=verbose,
            selmer_only=True,
            first_limit=first_limit,
            second_limit=second_limit,
            second_descent=second_descent,
        )

        upper = int(M.rank_bound())
        result["rank_upper_bound"] = upper
        result["rank_jump_upper"] = upper - int(generic_rank_reference)
        result["rank_jump_candidate"] = upper > int(generic_rank_reference)
        result["can_reach_target"] = upper >= int(target_rank)
        result["target_rejected_by_upper_bound"] = upper < int(target_rank)

    except Exception as exc:
        result["error"] = repr(exc)

    result["elapsed_seconds"] = time.time() - started
    return result


def run_mwrank_full(
    Emin,
    first_limit,
    second_limit,
    second_descent,
    verbose,
):
    """
    Full eclib/mwrank 2-descent including global Selmer-class resolution.

    Use this only for candidates whose Selmer-only screening upper bound
    exceeds the known generic-section lower bound.

    This records only lower bound, upper bound, and certainty. It deliberately
    avoids gens(), regulator(), Silverman bounds, and explicit saturation.
    """
    started = time.time()

    result = {
        "enabled": True,
        "mode": "full_global_descent",
        "first_limit": first_limit,
        "second_limit": second_limit,
        "second_descent": second_descent,
        "rank_lower_bound": None,
        "rank_upper_bound": None,
        "certain": None,
        "minimal_model_ainvs": [str(x) for x in Emin.a_invariants()],
        "error": None,
        "elapsed_seconds": None,
    }

    try:
        ainvs = [int(x) for x in Emin.a_invariants()]
        M = mwrank_EllipticCurve(ainvs, verbose=verbose)

        M.two_descent(
            verbose=verbose,
            selmer_only=False,
            first_limit=first_limit,
            second_limit=second_limit,
            second_descent=second_descent,
        )

        try:
            result["rank_lower_bound"] = int(M.rank())
        except Exception as exc:
            result["rank_error"] = repr(exc)

        try:
            result["rank_upper_bound"] = int(M.rank_bound())
        except Exception as exc:
            result["rank_bound_error"] = repr(exc)

        try:
            result["certain"] = bool(M.certain())
        except Exception as exc:
            result["certain_error"] = repr(exc)

    except Exception as exc:
        result["error"] = repr(exc)

    result["elapsed_seconds"] = time.time() - started
    return result


# ---------------------------------------------------------------------------
# Candidate analysis
# ---------------------------------------------------------------------------

def analyze_candidate(
    u,
    candidate,
    rel_tol,
    point_search_bound,
    run_pari,
    run_saturation,
    run_mwrank,
    run_mwrank_full_flag,
    generic_rank_reference,
    target_rank,
    mwrank_first_limit,
    mwrank_second_limit,
    mwrank_second_descent,
    mwrank_verbose,
):
    started = time.time()

    t = Q(candidate["t"])
    A = A_polynomial(u, t)
    B = B_product(u, t)

    E = EllipticCurve(
        QQ,
        [
            0,
            QQ(2*A.numerator) / A.denominator,
            0,
            QQ(B.numerator) / B.denominator,
            0,
        ],
    )

    Emin = None
    minimal_model = None
    minimal_model_error = None

    try:
        Emin = E.global_minimal_model()
        minimal_model = [str(x) for x in Emin.a_invariants()]
    except Exception as exc:
        minimal_model_error = repr(exc)

    torsion = None
    torsion_error = None

    try:
        torsion = str(E.torsion_subgroup())
    except Exception as exc:
        torsion_error = repr(exc)

    section_result = recover_generic_sections(E, u, t, A, B)
    generic_points = section_result.pop("points", [])

    generic_span = None

    if generic_points:
        try:
            generic_span = numerical_span(
                E,
                generic_points,
                rel_tol=rel_tol,
            )
        except Exception as exc:
            generic_span = {
                "estimated_span": None,
                "error": repr(exc),
            }

    point_search_rows = []
    point_search_error = None
    extra_points = []
    combined_span = None

    if point_search_bound > 0:
        try:
            found = E.point_search(
                bound=point_search_bound,
            )

            found = unique_points(found)
            extra_points = extra_points_only(
                generic_points,
                found,
            )

            point_search_rows = [
                {
                    "x": str(P[0]),
                    "y": str(P[1]),
                    "height": str(P.height()),
                }
                for P in extra_points
            ]

            if extra_points:
                combined = unique_points(
                    generic_points + extra_points
                )

                combined_span = numerical_span(
                    E,
                    combined,
                    rel_tol=rel_tol,
                )

        except Exception as exc:
            point_search_error = repr(exc)

    pari_rank = None
    pari_error = None
    pari_elapsed = None

    if run_pari:
        pari_started = time.time()

        try:
            pari_rank = int(
                E.rank(algorithm="pari")
            )
        except Exception as exc:
            pari_error = repr(exc)

        pari_elapsed = time.time() - pari_started

    mwrank_result = {
        "enabled": False,
    }

    if run_mwrank:
        if Emin is None:
            mwrank_result = {
                "enabled": True,
                "error": (
                    "Global minimal model unavailable: "
                    f"{minimal_model_error}"
                ),
            }
        else:
            mwrank_result = run_mwrank_screen(
                Emin=Emin,
                generic_rank_reference=generic_rank_reference,
                target_rank=target_rank,
                first_limit=mwrank_first_limit,
                second_limit=mwrank_second_limit,
                second_descent=mwrank_second_descent,
                verbose=mwrank_verbose,
            )

    mwrank_full_result = {
        "enabled": False,
    }

    if run_mwrank_full_flag:
        if Emin is None:
            mwrank_full_result = {
                "enabled": True,
                "error": (
                    "Global minimal model unavailable: "
                    f"{minimal_model_error}"
                ),
            }
        else:
            mwrank_full_result = run_mwrank_full(
                Emin=Emin,
                first_limit=mwrank_first_limit,
                second_limit=mwrank_second_limit,
                second_descent=mwrank_second_descent,
                verbose=mwrank_verbose,
            )

    saturation = None

    if run_saturation and generic_points:
        try:
            points_for_sat = unique_points(
                generic_points + extra_points
            )

            sat_basis, sat_index, _ = E.saturation(
                points_for_sat
            )

            saturation = {
                "index": str(sat_index),
                "basis_count": len(sat_basis),
                "basis": [
                    {
                        "x": str(P[0]),
                        "y": str(P[1]),
                    }
                    for P in sat_basis
                ],
            }

        except Exception as exc:
            saturation = {
                "error": repr(exc),
            }

    elapsed = time.time() - started

    return {
        "source_score_rank": candidate.get("_source_rank"),
        "t": frac_str(t),
        "mestre_nagao_score": candidate.get("score"),
        "fixed_score": candidate.get("fixed_score"),

        "curve": {
            "A": frac_str(A),
            "B": frac_str(B),
            "ainvs": [str(x) for x in E.a_invariants()],
            "discriminant": str(E.discriminant()),
            "j_invariant": str(E.j_invariant()),
            "minimal_model_ainvs": minimal_model,
            "minimal_model_error": minimal_model_error,
            "torsion": torsion,
            "torsion_error": torsion_error,
        },

        "generic_sections": {
            **section_result,
            "recovered_count": len(generic_points),
            "span": generic_span,
        },

        "point_search": {
            "bound": point_search_bound,
            "extra_point_count": len(extra_points),
            "extra_points": point_search_rows,
            "error": point_search_error,
        },

        "combined_span": combined_span,

        "pari": {
            "enabled": run_pari,
            "rank": pari_rank,
            "error": pari_error,
            "elapsed_seconds": pari_elapsed,
        },

        "mwrank": mwrank_result,
        "mwrank_full": mwrank_full_result,

        "saturation": saturation,

        "elapsed_seconds": elapsed,
    }


# ---------------------------------------------------------------------------
# Input/output
# ---------------------------------------------------------------------------

def load_scored_candidates(path):
    with open(path) as f:
        obj = json.load(f)

    candidates = obj.get("candidates", [])

    for i, row in enumerate(candidates, start=1):
        row["_source_rank"] = i

    return obj, candidates


def parse_args():
    ap = argparse.ArgumentParser()

    ap.add_argument(
        "scores_json",
        help="JSON file produced by score_specializations.py.",
    )

    ap.add_argument(
        "--u",
        default="2/5",
    )

    ap.add_argument(
        "--top",
        type=int,
        default=5,
    )

    ap.add_argument(
        "--exclude-t",
        default="11860/97527",
        help=(
            "Comma-separated t values to exclude. "
            "Use an empty string to exclude nothing."
        ),
    )

    ap.add_argument(
        "--point-search-bound",
        type=int,
        default=0,
        help=(
            "Sage point_search bound. Default 0 disables it. "
            "The rank-19 benchmark showed this search to be ineffective "
            "at small/moderate bounds for the present family."
        ),
    )

    ap.add_argument(
        "--rel-tol",
        type=float,
        default=1e-12,
    )

    ap.add_argument(
        "--pari",
        action="store_true",
        help=(
            "Attempt E.rank(algorithm='pari') and record the result as a "
            "diagnostic, not as a standalone rank certificate."
        ),
    )

    ap.add_argument(
        "--saturate",
        action="store_true",
        help="Attempt Sage saturation of known points.",
    )

    ap.add_argument(
        "--generic-rank-reference",
        type=int,
        default=9,
        help=(
            "Published generic rank used only as a reference when measuring "
            "possible specialization jumps. Default: 9 for this EK family."
        ),
    )

    ap.add_argument(
        "--target-rank",
        type=int,
        default=20,
        help=(
            "Comparison target rank. A Selmer upper bound below this value "
            "rules out reaching that target in this screen. Default: 20."
        ),
    )

    ap.add_argument(
        "--mwrank",
        action="store_true",
        help=(
            "Run direct eclib/mwrank 2-descent on the global minimal model and stop before post-descent basis/saturation work."
        ),
    )

    ap.add_argument(
        "--mwrank-full",
        action="store_true",
        help=(
            "Run full eclib/mwrank global descent. Use only for candidates "
            "that pass the Selmer-only screen."
        ),
    )

    ap.add_argument(
        "--mwrank-first-limit",
        type=int,
        default=20,
        help=(
            "mwrank two_descent first_limit. Sage/eclib default is 20."
        ),
    )

    ap.add_argument(
        "--mwrank-second-limit",
        type=int,
        default=8,
        help=(
            "mwrank two_descent second_limit. Sage/eclib default is 8."
        ),
    )

    ap.add_argument(
        "--mwrank-no-second-descent",
        action="store_true",
        help=(
            "Disable eclib's second descent. For curves with 2-torsion, "
            "Sage documentation recommends leaving second descent enabled "
            "unless examining Selmer details."
        ),
    )

    ap.add_argument(
        "--mwrank-verbose",
        action="store_true",
        help="Print eclib/mwrank descent progress.",
    )

    ap.add_argument(
        "--output",
        default="analyzed_specializations.json",
    )

    return ap.parse_args()


def main():
    args = parse_args()

    if args.top <= 0:
        raise SystemExit("--top must be positive")

    source, candidates = load_scored_candidates(
        args.scores_json
    )

    excluded = set()

    if args.exclude_t.strip():
        excluded = {
            frac_str(Q(x.strip()))
            for x in args.exclude_t.split(",")
            if x.strip()
        }

    selected = [
        row
        for row in candidates
        if frac_str(Q(row["t"])) not in excluded
    ][:args.top]

    print("High-scoring specialization analyzer")
    print(f"source = {args.scores_json}")
    print(f"u = {args.u}")
    print(f"selected candidates = {len(selected)}")
    print(f"point-search bound = {args.point_search_bound}")
    print(f"PARI diagnostic enabled = {args.pari}")
    print(f"mwrank early-exit screening enabled = {args.mwrank}")
    print(f"mwrank full descent enabled = {args.mwrank_full}")
    print(f"generic-rank reference = {args.generic_rank_reference}")
    print(f"comparison target rank = {args.target_rank}")
    print(f"saturation enabled = {args.saturate}")

    if args.mwrank:
        print(
            "mwrank Selmer-only two_descent: "
            f"first_limit={args.mwrank_first_limit} "
            f"second_limit={args.mwrank_second_limit} "
            f"second_descent={not args.mwrank_no_second_descent}"
        )

    print()

    u = Q(args.u)
    rows = []

    for i, candidate in enumerate(selected, start=1):
        print(
            f"[{i}/{len(selected)}] "
            f"source_rank={candidate.get('_source_rank')} "
            f"t={candidate['t']} "
            f"score={candidate.get('score')} ...",
            flush=True,
        )

        row = analyze_candidate(
            u=u,
            candidate=candidate,
            rel_tol=args.rel_tol,
            point_search_bound=args.point_search_bound,
            run_pari=args.pari,
            run_saturation=args.saturate,
            run_mwrank=args.mwrank,
            run_mwrank_full_flag=args.mwrank_full,
            generic_rank_reference=args.generic_rank_reference,
            target_rank=args.target_rank,
            mwrank_first_limit=args.mwrank_first_limit,
            mwrank_second_limit=args.mwrank_second_limit,
            mwrank_second_descent=not args.mwrank_no_second_descent,
            mwrank_verbose=args.mwrank_verbose,
        )

        rows.append(row)

        generic_span = None

        gs = row["generic_sections"].get("span")

        if isinstance(gs, dict):
            generic_span = gs.get("estimated_span")

        combined_span = None

        if isinstance(row["combined_span"], dict):
            combined_span = row["combined_span"].get(
                "estimated_span"
            )

        mw = row["mwrank"]
        mwf = row["mwrank_full"]

        print(
            f"    sections={row['generic_sections']['recovered_count']}/9 "
            f"generic_span≈{generic_span} "
            f"extra_points={row['point_search']['extra_point_count']} "
            f"combined_span≈{combined_span} "
            f"PARI_diag={row['pari']['rank']} "
            f"mwrank_upper={mw.get('rank_upper_bound')} "
            f"jump_upper={mw.get('rank_jump_upper')} "
            f"jump_candidate={mw.get('rank_jump_candidate')} "
            f"can_reach_target={mw.get('can_reach_target')} "
            f"full=[{mwf.get('rank_lower_bound')},"
            f"{mwf.get('rank_upper_bound')}] "
            f"full_certain={mwf.get('certain')} "
            f"time={row['elapsed_seconds']:.2f}s"
        )

    output_path = resolve_output_path(
        args.output,
        "data/analyzed",
    )

    payload = {
        "source_scores_file": str(args.scores_json),

        "provenance": {
            "primary_reference": {
                "authors": "Noam D. Elkies; Zev Klagsbrun",
                "title": (
                    "New Rank Records For Elliptic Curves "
                    "Having Rational Torsion"
                ),
                "year": 2020,
                "arxiv": "2003.00077",
                "doi": "10.2140/obs.2020.4.233",
            },
            "software_method": {
                "name": "eclib/mwrank via SageMath",
                "purpose": (
                    "Selmer-only 2-descent upper-bound screening with early exit"
                ),
            },
        },

        "analysis": {
            "u": frac_str(u),
            "top_requested": args.top,
            "excluded_t": sorted(excluded),
            "point_search_bound": args.point_search_bound,
            "rel_tol": args.rel_tol,
            "pari_enabled": args.pari,
            "saturation_enabled": args.saturate,
            "mwrank_enabled": args.mwrank,
            "mwrank_full_enabled": args.mwrank_full,
            "generic_rank_reference": args.generic_rank_reference,
            "target_rank": args.target_rank,
            "mwrank_first_limit": args.mwrank_first_limit,
            "mwrank_second_limit": args.mwrank_second_limit,
            "mwrank_second_descent": (
                not args.mwrank_no_second_descent
            ),
        },

        "candidates": rows,
    }

    with output_path.open("w") as f:
        json.dump(
            payload,
            f,
            indent=2,
        )

    print()
    print("Summary:")

    for row in rows:
        gs = row["generic_sections"].get("span")
        generic_span = (
            gs.get("estimated_span")
            if isinstance(gs, dict)
            else None
        )

        mw = row["mwrank"]
        mwf = row["mwrank_full"]

        print(
            f"t={row['t']} "
            f"score={row['mestre_nagao_score']} "
            f"sections={row['generic_sections']['recovered_count']}/9 "
            f"generic_span≈{generic_span} "
            f"mwrank_upper={mw.get('rank_upper_bound')} "
            f"jump_upper={mw.get('rank_jump_upper')} "
            f"jump_candidate={mw.get('rank_jump_candidate')} "
            f"can_reach_target={mw.get('can_reach_target')} "
            f"full_lower={mwf.get('rank_lower_bound')} "
            f"full_upper={mwf.get('rank_upper_bound')} "
            f"full_certain={mwf.get('certain')}"
        )

    print()
    print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()

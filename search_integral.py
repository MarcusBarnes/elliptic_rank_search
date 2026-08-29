#!/usr/bin/env python3
"""
Search for elliptic curves

    y^2 + x*y = x^3 + A*x + B

with many visible integral points constrained by

    x = 10U + 4
    y = 10V + 3.

Each (U,V) induces the auxiliary point

    (x, q(U,V))

where q(U,V) = A*x + B.

Thus any affine line through many auxiliary points gives one candidate curve.

This first implementation prioritizes transparency and reproducibility over scale.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from dataclasses import dataclass, asdict
from fractions import Fraction
from math import gcd
from typing import Dict, Iterable, List, Tuple
from pathlib import Path


@dataclass(frozen=True)
class AuxPoint:
    U: int
    V: int
    x: int
    y: int
    q: int


def q_value(U: int, V: int) -> int:
    return (
        100 * V * V
        + 100 * (U + 1) * V
        - 1000 * U**3
        - 1200 * U**2
        - 450 * U
        - 43
    )


def make_aux_points(umin: int, umax: int, vmin: int, vmax: int) -> List[AuxPoint]:
    pts: List[AuxPoint] = []
    for U in range(umin, umax + 1):
        x = 10 * U + 4
        for V in range(vmin, vmax + 1):
            y = 10 * V + 3
            pts.append(AuxPoint(U=U, V=V, x=x, y=y, q=q_value(U, V)))
    return pts


def normalize_line(p1: AuxPoint, p2: AuxPoint):
    """
    Return an exact line key (A,B), represented as Fractions, for q = A x + B.
    Points with equal x do not determine an affine line uniquely.
    """
    if p1.x == p2.x:
        return None
    A = Fraction(p2.q - p1.q, p2.x - p1.x)
    B = Fraction(p1.q) - A * p1.x
    return A, B


def discriminant_generalized(A: Fraction, B: Fraction) -> Fraction:
    """
    Discriminant for y^2 + x*y = x^3 + A*x + B,
    i.e. [a1,a2,a3,a4,a6] = [1,0,0,A,B].

    b2 = 1
    b4 = 2A
    b6 = 4B
    b8 = B - A^2

    Delta = -b2^2*b8 - 8*b4^3 - 27*b6^2 + 9*b2*b4*b6
    """
    b2 = Fraction(1)
    b4 = 2 * A
    b6 = 4 * B
    b8 = B - A * A
    return -(b2 * b2) * b8 - 8 * b4**3 - 27 * b6**2 + 9 * b2 * b4 * b6


def verify_point(A: Fraction, B: Fraction, x: int, y: int) -> bool:
    return Fraction(y*y + x*y) == Fraction(x**3) + A*x + B


def search(
    aux_points: List[AuxPoint],
    min_points: int,
    integer_coefficients_only: bool = True,
):
    """
    Enumerate lines determined by pairs of auxiliary points.

    We accumulate the set of point indices lying on each line.
    For correctness, after pair accumulation we rescan only candidate lines
    having enough pair support and recover every point on the line.

    This is O(N^2) and intended as a prototype.
    """
    line_support: Dict[Tuple[Fraction, Fraction], set] = defaultdict(set)
    n = len(aux_points)

    for i in range(n):
        p1 = aux_points[i]
        for j in range(i + 1, n):
            p2 = aux_points[j]
            if p1.x == p2.x:
                continue
            key = normalize_line(p1, p2)
            if key is None:
                continue
            A, B = key
            if integer_coefficients_only and (A.denominator != 1 or B.denominator != 1):
                continue
            line_support[key].add(i)
            line_support[key].add(j)

    results = []
    for (A, B), idxs in line_support.items():
        # pair support is only a lower bound; recover all actual points
        if len(idxs) < min_points:
            continue

        members = []
        seen_x = set()
        for p in aux_points:
            if Fraction(p.q) == A * p.x + B:
                # Rank comes from distinct elliptic-curve points. We score by
                # distinct x first to avoid counting both y-roots as excessive evidence.
                members.append(p)
                seen_x.add(p.x)

        if len(seen_x) < min_points:
            continue

        disc = discriminant_generalized(A, B)
        if disc == 0:
            continue

        assert all(verify_point(A, B, p.x, p.y) for p in members)

        results.append({
            "A": str(A),
            "B": str(B),
            "discriminant": str(disc),
            "distinct_x_count": len(seen_x),
            "visible_point_count": len(members),
            "points": [
                {"U": p.U, "V": p.V, "x": p.x, "y": p.y}
                for p in members
            ],
        })

    results.sort(
        key=lambda r: (
            r["distinct_x_count"],
            r["visible_point_count"],
            -len(r["A"]),
            -len(r["B"]),
        ),
        reverse=True,
    )
    return results


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--umin", type=int, default=-8)
    ap.add_argument("--umax", type=int, default=8)
    ap.add_argument("--vmin", type=int, default=-50)
    ap.add_argument("--vmax", type=int, default=50)
    ap.add_argument("--min-points", type=int, default=6)
    ap.add_argument("--top", type=int, default=20)
    ap.add_argument("--allow-rational-coefficients", action="store_true")
    ap.add_argument("--output", type=str)
    return ap.parse_args()
    
def resolve_output_path(filename, default_dir):
    path = Path(filename)

    # Bare filenames go into the standard project output directory.
    # Explicit paths are respected as supplied.
    if path.parent == Path("."):
        path = Path(default_dir) / path

    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def main():
    args = parse_args()
    aux = make_aux_points(args.umin, args.umax, args.vmin, args.vmax)
    print(f"Generated {len(aux)} auxiliary points")

    results = search(
        aux,
        min_points=args.min_points,
        integer_coefficients_only=not args.allow_rational_coefficients,
    )

    print(f"Found {len(results)} nonsingular candidate curves")
    for i, r in enumerate(results[: args.top]):
        print(
            f"[{i}] A={r['A']} B={r['B']} "
            f"distinct_x={r['distinct_x_count']} visible_points={r['visible_point_count']}"
        )

    if args.output:
        payload = {
            "search": {
                "umin": args.umin,
                "umax": args.umax,
                "vmin": args.vmin,
                "vmax": args.vmax,
                "min_points": args.min_points,
                "integer_coefficients_only": not args.allow_rational_coefficients,
            },
            "candidates": results[: args.top],
        }
        output_path = resolve_output_path(
            args.output,
            "data/candidates",
        )

        with output_path.open("w") as f:
            json.dump(payload, f, indent=2)

        print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()

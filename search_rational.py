#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction
from typing import Dict, Iterable, List, Tuple
from pathlib import Path


@dataclass(frozen=True)
class AuxPoint:
    U: int
    V: int
    Z: int
    X: int
    Y: int
    x: Fraction
    y: Fraction
    q: Fraction

def resolve_output_path(filename, default_dir):
    path = Path(filename)

    if path.parent == Path("."):
        path = Path(default_dir) / path

    path.parent.mkdir(parents=True, exist_ok=True)
    return path

def frac_to_str(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def parse_z_values(text: str) -> List[int]:
    vals = []
    for item in text.split(","):
        item = item.strip()
        if not item:
            continue
        z = int(item)
        if z <= 0:
            raise argparse.ArgumentTypeError("All Z values must be positive integers.")
        vals.append(z)
    if not vals:
        raise argparse.ArgumentTypeError("At least one Z value is required.")
    return list(dict.fromkeys(vals))


def make_aux_point(U: int, V: int, Z: int) -> AuxPoint:
    Z2 = Z * Z
    Z3 = Z2 * Z
    Z6 = Z3 * Z3
    X = 4 * Z2 + 10 * U
    Y = 3 * Z3 + 10 * V
    x = Fraction(X, Z2)
    y = Fraction(Y, Z3)
    q = Fraction(Y * Y + X * Y * Z - X**3, Z6)
    assert q == y * y + x * y - x**3
    return AuxPoint(U, V, Z, X, Y, x, y, q)


def make_aux_points(
    umin: int,
    umax: int,
    vmin: int,
    vmax: int,
    z_values: Iterable[int],
) -> List[AuxPoint]:
    unique: Dict[Tuple[Fraction, Fraction], AuxPoint] = {}
    for Z in z_values:
        for U in range(umin, umax + 1):
            for V in range(vmin, vmax + 1):
                p = make_aux_point(U, V, Z)
                unique.setdefault((p.x, p.y), p)
    return list(unique.values())


def normalize_line(p1: AuxPoint, p2: AuxPoint):
    if p1.x == p2.x:
        return None
    A = (p2.q - p1.q) / (p2.x - p1.x)
    B = p1.q - A * p1.x
    return A, B


def discriminant_generalized(A: Fraction, B: Fraction) -> Fraction:
    b2 = Fraction(1)
    b4 = 2 * A
    b6 = 4 * B
    b8 = B - A * A
    return -(b2 * b2) * b8 - 8 * b4**3 - 27 * b6**2 + 9 * b2 * b4 * b6


def verify_point(A: Fraction, B: Fraction, p: AuxPoint) -> bool:
    return p.y * p.y + p.x * p.y == p.x**3 + A * p.x + B


def search(
    aux_points: List[AuxPoint],
    min_points: int,
    integer_coefficients_only: bool = False,
):
    line_support = defaultdict(set)
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
            if integer_coefficients_only and (
                A.denominator != 1 or B.denominator != 1
            ):
                continue
            line_support[key].add(i)
            line_support[key].add(j)

    results = []
    for (A, B), idxs in line_support.items():
        if len(idxs) < min_points:
            continue

        members = []
        seen_x = set()
        for p in aux_points:
            if p.q == A * p.x + B:
                members.append(p)
                seen_x.add(p.x)

        if len(seen_x) < min_points:
            continue

        disc = discriminant_generalized(A, B)
        if disc == 0:
            continue

        assert all(verify_point(A, B, p) for p in members)

        results.append({
            "A": frac_to_str(A),
            "B": frac_to_str(B),
            "discriminant": frac_to_str(disc),
            "distinct_x_count": len(seen_x),
            "visible_point_count": len(members),
            "points": [
                {
                    "U": p.U,
                    "V": p.V,
                    "Z": p.Z,
                    "X": p.X,
                    "Y": p.Y,
                    "x": frac_to_str(p.x),
                    "y": frac_to_str(p.y),
                }
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
    ap.add_argument("--umin", type=int, default=-6)
    ap.add_argument("--umax", type=int, default=6)
    ap.add_argument("--vmin", type=int, default=-30)
    ap.add_argument("--vmax", type=int, default=30)
    ap.add_argument(
        "--z-values",
        type=parse_z_values,
        default=parse_z_values("1,2,3"),
        help="Comma-separated positive Z values, e.g. 1,2,3,5.",
    )
    ap.add_argument("--min-points", type=int, default=6)
    ap.add_argument("--top", type=int, default=100)
    ap.add_argument("--integer-coefficients-only", action="store_true")
    ap.add_argument("--output", type=str)
    return ap.parse_args()


def main():
    args = parse_args()

    aux = make_aux_points(
        args.umin,
        args.umax,
        args.vmin,
        args.vmax,
        args.z_values,
    )

    print(
        f"Generated {len(aux)} unique rational auxiliary points "
        f"for Z={args.z_values}"
    )

    results = search(
        aux,
        min_points=args.min_points,
        integer_coefficients_only=args.integer_coefficients_only,
    )

    print(f"Found {len(results)} nonsingular candidate curves")

    for i, r in enumerate(results[: args.top]):
        zs = sorted({p["Z"] for p in r["points"]})
        print(
            f"[{i}] A={r['A']} B={r['B']} "
            f"distinct_x={r['distinct_x_count']} "
            f"visible_points={r['visible_point_count']} "
            f"Z={zs}"
        )

    if args.output:
        payload = {
            "search": {
                "type": "rational_residue_class",
                "umin": args.umin,
                "umax": args.umax,
                "vmin": args.vmin,
                "vmax": args.vmax,
                "z_values": args.z_values,
                "min_points": args.min_points,
                "integer_coefficients_only": args.integer_coefficients_only,
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

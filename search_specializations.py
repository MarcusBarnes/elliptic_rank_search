#!/usr/bin/env sage
from __future__ import annotations

import argparse
import json
from fractions import Fraction
from math import isqrt
from pathlib import Path

def resolve_output_path(filename, default_dir):
    path = Path(filename)

    if path.parent == Path("."):
        path = Path(default_dir) / path

    path.parent.mkdir(parents=True, exist_ok=True)
    return path

def Q(s) -> Fraction:
    if isinstance(s, Fraction):
        return s
    return Fraction(str(s))


def frac_str(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def rational_square_root(x: Fraction):
    if x < 0:
        return None
    ra = isqrt(x.numerator)
    rb = isqrt(x.denominator)
    if ra * ra == x.numerator and rb * rb == x.denominator:
        return Fraction(ra, rb)
    return None


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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--u", required=True)
    ap.add_argument("--t", required=True)
    ap.add_argument("--output")
    ap.add_argument(
        "--saturate",
        action="store_true",
        help="Run Sage saturation on the recovered sections (expensive).",
    )
    args = ap.parse_args()

    from sage.all import QQ, RR, EllipticCurve

    u = Q(args.u)
    t = Q(args.t)

    A = A_polynomial(u, t)
    bfs = B_factors(u, t)
    B = Fraction(1)
    for f in bfs:
        B *= f

    E = EllipticCurve(QQ, [0, QQ(2*A.numerator)/A.denominator, 0,
                           QQ(B.numerator)/B.denominator, 0])

    square_term = Fraction(5) - u*u
    sqrt_s = rational_square_root(square_term)

    result = {
        "u": frac_str(u),
        "t": frac_str(t),
        "A": frac_str(A),
        "B": frac_str(B),
        "nonsingular": E.discriminant() != 0,
        "torsion": str(E.torsion_subgroup()),
        "minimal_model_ainvs": [str(x) for x in E.global_minimal_model().a_invariants()],
        "m_candidates": [],
        "sections": None,
    }

    print("Elkies--Klagsbrun specialization")
    print(f"u = {u}")
    print(f"t = {t}")
    print()
    print(f"5-u^2 = {square_term}")
    print(f"rational square? {sqrt_s is not None}")
    print(f"sqrt(5-u^2) = {sqrt_s}")
    print()

    if sqrt_s is None:
        print("Cannot construct ninth section: 5-u^2 is not a rational square.")
    else:
        m_pairs = compatible_m_values(u, sqrt_s)
        result["m_candidates"] = [
            {"m": frac_str(m), "signed_sqrt": frac_str(s)}
            for m, s in m_pairs
        ]

        print("Compatible m values:")
        for m, s in m_pairs:
            print(f"  m={m} signed_sqrt={s}")

        if not m_pairs:
            print("No compatible rational m found.")
        else:
            m_pairs.sort(key=lambda pair: pair[1] != sqrt_s)
            m = m_pairs[0][0]
            xs = published_section_xs(u, t, m)

            points = []
            rows = []
            failures = []

            for i, xf in enumerate(xs, start=1):
                # Recover y directly using exact rational arithmetic.
                # This avoids Sage's generic lift_x(), which can become
                # extremely expensive for the very large rational numbers
                # arising in record-rank specializations.
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

            sec = {
                "selected_m": frac_str(m),
                "point_count": len(points),
                "failed_sections": failures,
                "points": rows,
                "estimated_span": None,
                "height_determinant": None,
                "regulator_of_points": None,
                "saturation_index": None,
                "saturation_basis": None,
                "saturation_error": None,
            }

            if points:
                H = E.height_pairing_matrix(points)
                eigs = sorted([abs(RR(v)) for v in H.eigenvalues()], reverse=True)
                scale = max(eigs) if eigs else RR(0)
                tol = scale * RR("1e-12") if scale else RR("1e-12")
                sec["estimated_span"] = int(sum(1 for v in eigs if v > tol))
                sec["height_determinant"] = str(H.det())
                try:
                    sec["regulator_of_points"] = str(E.regulator_of_points(points))
                except Exception:
                    pass

                if args.saturate:
                    try:
                        sat_basis, sat_index, _ = E.saturation(points)
                        sec["saturation_index"] = str(sat_index)
                        sec["saturation_basis"] = [
                            [str(P[0]), str(P[1])] for P in sat_basis
                        ]
                    except Exception as exc:
                        sec["saturation_error"] = repr(exc)
                else:
                    sec["saturation_error"] = "Skipped (use --saturate to enable)"

            result["sections"] = sec

            print()
            print(f"Published sections recovered: {len(points)}/9")
            print(f"selected m = {m}")
            print(f"estimated section span ≈ {sec['estimated_span']}")
            print(f"height determinant = {sec['height_determinant']}")
            print(f"regulator of points = {sec['regulator_of_points']}")
            print(f"saturation index = {sec['saturation_index']}")

            if failures:
                print("Failed sections:")
                for f in failures:
                    print(f"  P{f['section']}: x={f['x']} {f['error']}")

            print()
            print("Specialized section points:")
            for r in rows:
                print(f"  P{r['section']}: x={r['x']} y={r['y']}")

    print()
    print(f"torsion = {result['torsion']}")
    print(f"minimal model ainvs = {result['minimal_model_ainvs']}")

    if args.output:
        output_path = resolve_output_path(
            args.output,
            "data/benchmarks",
        )
    
        with output_path.open("w") as f:
            json.dump(result, f, indent=2)
    
        print()
        print(f"Wrote {output_path}")

if __name__ == "__main__":
    main()

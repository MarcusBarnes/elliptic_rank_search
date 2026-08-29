#!/usr/bin/env sage
from __future__ import annotations

import argparse
import csv
import json
import time
from pathlib import Path

from sage.all import QQ, EllipticCurve


def resolve_output_path(filename, default_dir):
    path = Path(filename)

    if path.parent == Path("."):
        path = Path(default_dir) / path

    path.parent.mkdir(parents=True, exist_ok=True)
    return path

def q(value):
    return QQ(value)


def load_candidates(path):
    with open(path) as f:
        obj = json.load(f)
    if "candidates" in obj:
        return obj["candidates"], obj.get("search", {})
    return [obj], {}


def one_point_per_x(points):
    chosen = {}
    for p in points:
        x = q(p["x"])
        if x not in chosen:
            chosen[x] = p
    return list(chosen.values())


def denominator_usage(points):
    """Summarize Z usage. Integral-search points without Z are treated as Z=1."""
    z_values = set()
    nonintegral_point_count = 0
    nonintegral_x = set()

    for p in points:
        z = int(p.get("Z", 1))
        z_values.add(z)
        if z > 1:
            nonintegral_point_count += 1
            nonintegral_x.add(q(p["x"]))

    if not z_values:
        z_values.add(1)

    return sorted(z_values), nonintegral_point_count, len(nonintegral_x)


def height_span_estimate(E, points, rel_tol=1e-10):
    if not points:
        return 0, [], None

    H = E.height_pairing_matrix(points)
    eigs = sorted([abs(x) for x in H.eigenvalues()], reverse=True)

    if not eigs:
        return 0, [], H.det()

    scale = max(eigs)
    tol = rel_tol if scale == 0 else scale * rel_tol
    estimated_span = sum(1 for x in eigs if x > tol)
    return estimated_span, eigs, H.det()


def screen_candidate(candidate, index, rel_tol=1e-10):
    started = time.time()

    A = q(candidate["A"])
    B = q(candidate["B"])
    E = EllipticCurve(QQ, [1, 0, 0, A, B])

    raw_points = candidate.get("points", [])
    z_values_used, nonintegral_point_count, nonintegral_distinct_x_count = denominator_usage(raw_points)

    representative_dicts = one_point_per_x(raw_points)
    reps = [E(q(p["x"]), q(p["y"])) for p in representative_dicts]

    estimated_span = None
    eigs = []
    height_det = None
    height_error = None

    try:
        estimated_span, eigs, height_det = height_span_estimate(E, reps, rel_tol=rel_tol)
    except Exception as exc:
        height_error = repr(exc)

    torsion_order = None
    torsion_description = None
    try:
        T = E.torsion_subgroup()
        torsion_order = int(T.order())
        torsion_description = str(T)
    except Exception as exc:
        torsion_description = f"ERROR: {exc!r}"

    return {
        "index": index,
        "A": str(A),
        "B": str(B),
        "discriminant": str(E.discriminant()),
        "j_invariant": str(E.j_invariant()),
        "raw_visible_points": len(raw_points),
        "distinct_x_points": len(reps),
        "z_values_used": z_values_used,
        "nonintegral_point_count": nonintegral_point_count,
        "nonintegral_distinct_x_count": nonintegral_distinct_x_count,
        "estimated_span": estimated_span,
        "height_eigenvalues": [str(x) for x in eigs],
        "height_determinant": None if height_det is None else str(height_det),
        "height_error": height_error,
        "pari_rank": None,
        "pari_error": None,
        "analytic_rank": None,
        "analytic_error": None,
        "torsion_order": torsion_order,
        "torsion_description": torsion_description,
        "screen_elapsed_seconds": time.time() - started,
        "rank_elapsed_seconds": None,
    }


def compute_ranks(row, analytic=False):
    started = time.time()
    A = q(row["A"])
    B = q(row["B"])
    E = EllipticCurve(QQ, [1, 0, 0, A, B])

    try:
        row["pari_rank"] = int(E.rank(algorithm="pari"))
    except Exception as exc:
        row["pari_error"] = repr(exc)

    if analytic:
        try:
            row["analytic_rank"] = int(E.analytic_rank())
        except Exception as exc:
            row["analytic_error"] = repr(exc)

    row["rank_elapsed_seconds"] = time.time() - started


def screening_sort_key(row):
    span = row["estimated_span"] if row["estimated_span"] is not None else -1
    return (span, row["distinct_x_points"], row["nonintegral_distinct_x_count"])


def final_sort_key(row):
    pari = row["pari_rank"] if row["pari_rank"] is not None else -1
    span = row["estimated_span"] if row["estimated_span"] is not None else -1
    return (pari, span, row["distinct_x_points"], row["nonintegral_distinct_x_count"])


def write_csv(path, rows):
    fields = [
        "index", "A", "B",
        "raw_visible_points", "distinct_x_points",
        "z_values_used", "nonintegral_point_count", "nonintegral_distinct_x_count",
        "estimated_span", "pari_rank", "analytic_rank", "torsion_order",
        "discriminant", "height_determinant",
        "screen_elapsed_seconds", "rank_elapsed_seconds",
    ]

    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            csv_row = {k: row.get(k) for k in fields}
            csv_row["z_values_used"] = ",".join(str(z) for z in row.get("z_values_used", []))
            writer.writerow(csv_row)


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("json_file")
    ap.add_argument(
        "--output-json",
        default="data/analyzed/analyzed_candidates.json",
    )
    
    ap.add_argument(
        "--output-csv",
        default="data/analyzed/analyzed_candidates.csv",
    )
    ap.add_argument("--analytic", action="store_true")
    ap.add_argument("--rel-tol", type=float, default=1e-10)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument(
        "--pari-top-k",
        type=int,
        default=100,
        help="Run PARI only on the top K screened candidates; use 0 for all.",
    )
    return ap.parse_args()


def main():
    args = parse_args()
    candidates, search_metadata = load_candidates(args.json_file)

    if args.limit is not None:
        candidates = candidates[:args.limit]

    total = len(candidates)
    print(f"Stage 1: screening {total} candidate(s)")
    print(f"Numerical span relative tolerance: {args.rel_tol:g}")

    rows = []
    for i, candidate in enumerate(candidates):
        print(
            f"[{i + 1}/{total}] A={candidate['A']} B={candidate['B']} ...",
            end=" ", flush=True,
        )
        row = screen_candidate(candidate, index=i, rel_tol=args.rel_tol)
        rows.append(row)
        print(
            f"distinct_x={row['distinct_x_points']} "
            f"span≈{row['estimated_span']} "
            f"Z={row['z_values_used']} "
            f"nonint_x={row['nonintegral_distinct_x_count']} "
            f"time={row['screen_elapsed_seconds']:.2f}s"
        )

    rows.sort(key=screening_sort_key, reverse=True)
    shortlist_size = len(rows) if args.pari_top_k == 0 else min(args.pari_top_k, len(rows))

    print()
    print(f"Stage 2: running PARI rank on {shortlist_size} candidate(s)")

    for i, row in enumerate(rows[:shortlist_size]):
        print(
            f"[{i + 1}/{shortlist_size}] index={row['index']} "
            f"A={row['A']} B={row['B']} span≈{row['estimated_span']} "
            f"Z={row['z_values_used']} ...",
            end=" ", flush=True,
        )
        compute_ranks(row, analytic=args.analytic)
        print(
            f"PARI={row['pari_rank']} analytic={row['analytic_rank']} "
            f"time={row['rank_elapsed_seconds']:.2f}s"
        )

    rows.sort(key=final_sort_key, reverse=True)

    payload = {
        "source_file": str(args.json_file),
        "source_search": search_metadata,
        "analysis": {
            "rel_tol": args.rel_tol,
            "analytic_enabled": args.analytic,
            "candidate_count": len(rows),
            "pari_top_k": shortlist_size,
        },
        "candidates": rows,
    }
    
    
    output_json = resolve_output_path(
        args.output_json,
        "data/analyzed",
    )
    
    output_csv = resolve_output_path(
        args.output_csv,
        "data/analyzed",
    )

    with output_json.open("w") as f:
        json.dump(payload, f, indent=2)

    write_csv(
        output_csv,
        rows,
    )

    print()
    print("Top candidates:")
    for row in rows[:20]:
        print(
            f"index={row['index']} A={row['A']} B={row['B']} "
            f"distinct_x={row['distinct_x_points']} "
            f"span≈{row['estimated_span']} "
            f"PARI={row['pari_rank']} analytic={row['analytic_rank']} "
            f"Z={row['z_values_used']} "
            f"nonint_points={row['nonintegral_point_count']} "
            f"nonint_x={row['nonintegral_distinct_x_count']}"
        )

    print()
    print(f"Wrote {output_json}")
    print(f"Wrote {output_csv}")


if __name__ == "__main__":
    main()

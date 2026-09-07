#!/usr/bin/env sage
"""
Reproduce Fermigier's published 1997 rank-at-least-22 specialization.

This is a literature-reproduction artifact, not a new rank claim and not a
specialization-search program.  It verifies:

1. the normalized split sextic used for the selected Mestre--Fermigier family;
2. the parameter convention t_paper = 19754/39 and symmetric shift T=2*t_paper;
3. exact square completion q6(x-T)q6(x+T)=g(x)^2-r(x), deg(r)=4;
4. the twelve forced rational x-coordinates on y^2=r(x);
5. QQ-isomorphism of the quartic Jacobian with Fermigier's published E22 model;
6. Fermigier's published Mestre--Nagao score checkpoints.

Primary source:
Stéfane Fermigier, "Une courbe elliptique définie sur Q de rang >= 22",
Acta Arithmetica 82(4) (1997), 359--363.
DOI: 10.4064/aa-82-4-359-363.

The paper's rank >= 22 result is cited, not re-certified by this script.
"""

from sage.all import QQ, ZZ, RR, PolynomialRing, EllipticCurve, Jacobian
from pathlib import Path
from datetime import datetime
import json
import time

SELECTED = [QQ(0), QQ(55), QQ(314), QQ(378), QQ(1007), QQ(1036)]
TPAPER = QQ(19754) / QQ(39)
T = 2 * TPAPER

PUBLISHED_E22 = EllipticCurve(
    QQ,
    [
        1, 0, 1,
        ZZ(-940299517776391362903023121165864),
        ZZ(10707363070719743033425295515449274534651125011362),
    ],
)

PUBLISHED_SCORES = {
    50: 29.49,
    100: 44.12,
    200: 57.54,
    400: 81.51,
    1000: 105.17,
    2000: 122.76,
}


def polynomial_part_square_root(q, x):
    if q.degree() != 12 or q[12] != 1:
        raise RuntimeError("expected a monic degree-12 polynomial")
    g = x**6
    for d in range(5, -1, -1):
        degree = 6 + d
        current = (g * g)[degree]
        target = q[degree]
        g += ((target - current) / 2) * x**d
    r = g * g - q
    return g, r


def first_n_primes(n):
    out = []
    p = ZZ(1)
    for _ in range(n):
        p = p.next_prime()
        out.append(p)
    return out


def historical_score(E, M, primes):
    """Fermigier's printed convention: first M primes, omitting p=2."""
    total = RR(0)
    for p in primes[:M]:
        if p == 2:
            continue
        ap = ZZ(E.ap(p))
        Np = ZZ(p + 1 - ap)
        total += RR(2 - ap) / RR(Np) * RR(p).log()
    return float(total)


def main():
    started = datetime.now().astimezone().isoformat()
    wall0 = time.perf_counter()

    Rx = PolynomialRing(QQ, "x")
    x = Rx.gen()
    q6 = Rx(1)
    for a in SELECTED:
        q6 *= x - a

    q12 = q6(x - T) * q6(x + T)
    g, r = polynomial_part_square_root(q12, x)

    identity_ok = (q12 == g * g - r)
    quartic_ok = (r.degree() == 4)

    forced = []
    forced_ok = True
    for a in SELECTED:
        for sign in (-1, 1):
            xx = a + sign * T
            yy = g(xx)
            ok = (yy * yy == r(xx))
            forced_ok = forced_ok and ok
            forced.append({"x": str(xx), "y": str(yy), "verified": bool(ok)})

    Rxy = PolynomialRing(QQ, names=("X", "Y"))
    X, Y = Rxy.gens()
    rXY = sum(QQ(r[i]) * X**i for i in range(5))
    E = Jacobian(Y**2 - rXY)

    j_match = (E.j_invariant() == PUBLISHED_E22.j_invariant())
    isomorphic = E.is_isomorphic(PUBLISHED_E22)

    primes = first_n_primes(max(PUBLISHED_SCORES))
    score_rows = []
    scores_ok = True
    for M, published in PUBLISHED_SCORES.items():
        computed = historical_score(E, M, primes)
        # Historical table is printed to two decimal places.
        match = (round(computed, 2) == round(published, 2))
        scores_ok = scores_ok and match
        score_rows.append({
            "M": int(M),
            "computed": float(computed),
            "published": float(published),
            "match_2dp": bool(match),
        })

    passed = identity_ok and quartic_ok and forced_ok and j_match and isomorphic and scores_ok

    payload = {
        "source": {
            "author": "Stéfane Fermigier",
            "title": "Une courbe elliptique définie sur Q de rang >= 22",
            "journal": "Acta Arithmetica",
            "volume": "82",
            "issue": "4",
            "year": int(1997),
            "pages": "359-363",
            "doi": "10.4064/aa-82-4-359-363",
        },
        "proof_status": "published-result reproduction; rank >= 22 is not independently re-certified here",
        "selected_roots": [str(a) for a in SELECTED],
        "t_paper": str(TPAPER),
        "symmetric_shift_T": str(T),
        "identity_verified": bool(identity_ok),
        "quartic_degree": int(r.degree()),
        "forced_points_verified": bool(forced_ok),
        "forced_points": forced,
        "jacobian_ainvs": [str(a) for a in E.ainvs()],
        "published_E22_ainvs": [str(a) for a in PUBLISHED_E22.ainvs()],
        "j_invariant_match": bool(j_match),
        "QQ_isomorphic_to_published_E22": bool(isomorphic),
        "historical_score_convention": "first M primes, p=2 omitted",
        "historical_scores": score_rows,
        "all_checks_passed": bool(passed),
        "started_at": started,
        "wall_seconds": float(time.perf_counter() - wall0),
    }

    here = Path(__file__).resolve().parent
    out = here / "reproduction_result.json"
    out.write_text(json.dumps(payload, indent=2) + "\n")

    print("Fermigier 1997 rank >= 22 literature reproduction")
    print("=" * 72)
    print("t_paper:", TPAPER)
    print("symmetric shift T:", T)
    print("quartic degree:", r.degree())
    print("exact identity:", identity_ok)
    print("12 forced points verified:", forced_ok)
    print("j-invariant match:", j_match)
    print("QQ-isomorphic to published E22:", isomorphic)
    print()
    print("Historical score checkpoints:")
    for row in score_rows:
        print(
            "  M={:4d}: computed={:10.4f} published={:6.2f} match2dp={}".format(
                row["M"], row["computed"], row["published"], row["match_2dp"]
            )
        )
    print()
    print("ALL CHECKS PASSED:", passed)
    print("Result JSON:", out)
    print("wall seconds: {:.3f}".format(payload["wall_seconds"]))

    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

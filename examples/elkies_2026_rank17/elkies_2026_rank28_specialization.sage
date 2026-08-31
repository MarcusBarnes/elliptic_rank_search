#!/usr/bin/env sage
"""
Public-safe reproduction of Elkies's published rank-28 specialization.

Sources
-------
1. Noam D. Elkies (2026), arXiv:2608.25406v1, Section 2.3:
     the rank-at-least-28 fiber occurs at t = -9529/5471.
2. Klagsbrun--Sherman--Weigandt (2016), arXiv:1606.07178:
     published generalized Weierstrass model E_28.
3. Published reproductions of Elkies's 28 independent rational points,
   e.g. Dujella/Yokoyama lecture material derived from Elkies's 2006
   NMBRTHRY announcement.

Scope
-----
This script only reproduces published data:
  * specializes the published rank-17 K3 fibration at t = -9529/5471;
  * checks exact Q-isomorphism to the published E_28 model;
  * checks specialization of the 17 published generic sections;
  * checks the 28 published rational points on E_28;
  * computes their Neron--Tate height-pairing matrix and numerical regulator
    as a reproducibility check of independence.

It contains no unpublished specialization search, ranking heuristic,
candidate generation, or Track A/Track B extension.
"""

from pathlib import Path
from sage.all import QQ, ZZ, EllipticCurve, matrix, load

HERE = Path(__file__).resolve().parent

# Load and re-run the generic-rank-17 reproduction.  This defines
# t, A, B, xs, ys after its own exact checks have passed.
load(str(HERE / "elkies_2026_rank17.sage"))

# Elkies 2026, Section 2.3.
t28 = -QQ(9529) / QQ(5471)

A28_short = QQ(A(t28))
B28_short = QQ(B(t28))

# Fiber in the short Weierstrass model inherited from the K3 fibration:
#     y^2 = x^3 + A(t28) x + B(t28).
E_fiber = EllipticCurve(QQ, [0, 0, 0, A28_short, B28_short])

# Published Elkies rank-28 model, KSW 2016, equation (1):
# y^2 + x*y + y = x^3 - x^2 + a4*x + a6.
a4 = -ZZ(20067762415575526585033208209338542750930230312178956502)
a6 = ZZ(34481611795030556467032985690390720374855944359319180361266008296291939448732243429)
E28 = EllipticCurve(QQ, [1, -1, 1, a4, a6])

# Exact fiber identification.
assert E_fiber.j_invariant() == E28.j_invariant()
assert E_fiber.is_isomorphic(E28)

# The 17 generic sections specialize to rational points on the fiber.
specialized_generic_points = [
    E_fiber(QQ(x(t28)), QQ(y(t28)))
    for x, y in zip(xs, ys)
]
assert len(specialized_generic_points) == 17

# Published 28 independent points on E_28.
# These coordinates are reproduced in multiple public sources based on
# Elkies's 3 May 2006 NMBRTHRY announcement.
published_points_xy = [
    (-2124150091254381073292137463, 259854492051899599030515511070780628911531),
    ( 2334509866034701756884754537,  18872004195494469180868316552803627931531),
    (-1671736054062369063879038663, 251709377261144287808506947241319126049131),
    ( 2139130260139156666492982137,  36639509171439729202421459692941297527531),
    ( 1534706764467120723885477337,  85429585346017694289021032862781072799531),
    (-2731079487875677033341575063, 262521815484332191641284072623902143387531),
    ( 2775726266844571649705458537,  12845755474014060248869487699082640369931),
    ( 1494385729327188957541833817,  88486605527733405986116494514049233411451),
    ( 1868438228620887358509065257,  59237403214437708712725140393059358589131),
    ( 2008945108825743774866542537,  47690677880125552882151750781541424711531),
    ( 2348360540918025169651632937,  17492930006200557857340332476448804363531),
    (-1472084007090481174470008663, 246643450653503714199947441549759798469131),
    ( 2924128607708061213363288937,  28350264431488878501488356474767375899531),
    ( 5374993891066061893293934537, 286188908427263386451175031916479893731531),
    ( 1709690768233354523334008557,  71898834974686089466159700529215980921631),
    ( 2450954011353593144072595187,   4445228173532634357049262550610714736531),
    ( 2969254709273559167464674937,  32766893075366270801333682543160469687531),
    ( 2711914934941692601332882937,   2068436612778381698650413981506590613531),
    (20078586077996854528778328937, 2779608541137806604656051725624624030091531),
    ( 2158082450240734774317810697,  34994373401964026809969662241800901254731),
    ( 2004645458247059022403224937,  48049329780704645522439866999888475467531),
    ( 2975749450947996264947091337,  33398989826075322320208934410104857869131),
    (-2102490467686285150147347863, 259576391459875789571677393171687203227531),
    (  311583179915063034902194537, 168104385229980603540109472915660153473931),
    ( 2773931008341865231443771817,  12632162834649921002414116273769275813451),
    ( 2156581188143768409363461387,  35125092964022908897004150516375178087331),
    ( 3866330499872412508815659137, 121197755655944226293036926715025847322531),
    ( 2230868289773576023778678737,  28558760030597485663387020600768640028531),
]

published_points = [E28(ZZ(x), ZZ(y)) for x, y in published_points_xy]
assert len(published_points) == 28

# Numerical Neron--Tate pairing check.  A nonzero, positive determinant is
# strong reproducibility evidence for the published independence claim;
# the publication itself is the source for the mathematical claim.
H = E28.height_pairing_matrix(published_points, precision=128)
reg = H.det()

print()
print("Elkies published rank-28 specialization reproduction")
print("  t =", t28)
print("  specialized fiber j-invariant matches published E_28")
print("  specialized fiber is Q-isomorphic to published E_28")
print("  17 generic sections specialize to rational points")
print("  28 published rational points verified on E_28")
print("  height-pairing matrix size =", H.nrows(), "x", H.ncols())
print("  numerical regulator =", reg)
print("  regulator positive =", bool(reg > 0))

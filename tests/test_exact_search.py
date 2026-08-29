import argparse
import json
import unittest
from fractions import Fraction
from pathlib import Path

import search_integral
import search_rational


class IntegralSearchTests(unittest.TestCase):
    def test_q_value_known_values(self):
        self.assertEqual(search_integral.q_value(0, 0), -43)
        self.assertEqual(search_integral.q_value(1, -2), -2693)

    def test_example_curve_discriminant_and_points(self):
        path = Path("example_curve.json")
        with path.open() as f:
            curve = json.load(f)

        A = Fraction(curve["A"])
        B = Fraction(curve["B"])

        self.assertEqual(
            search_integral.discriminant_generalized(A, B),
            Fraction(curve["discriminant"]),
        )

        for point in curve["points"]:
            self.assertTrue(
                search_integral.verify_point(A, B, point["x"], point["y"])
            )

    def test_small_integral_search_regression(self):
        aux = search_integral.make_aux_points(-1, 1, -3, 3)
        results = search_integral.search(aux, min_points=3)

        self.assertGreaterEqual(len(results), 1)
        top = results[0]
        self.assertEqual(top["A"], "-115")
        self.assertEqual(top["B"], "417")
        self.assertEqual(top["discriminant"], "18776000")
        self.assertEqual(top["distinct_x_count"], 3)
        self.assertEqual(top["visible_point_count"], 5)

        A = Fraction(top["A"])
        B = Fraction(top["B"])
        for point in top["points"]:
            self.assertTrue(
                search_integral.verify_point(A, B, point["x"], point["y"])
            )


class RationalSearchTests(unittest.TestCase):
    def test_z_one_agrees_with_integral_auxiliary_formula(self):
        p = search_rational.make_aux_point(2, -3, 1)

        self.assertEqual(p.x, Fraction(24))
        self.assertEqual(p.y, Fraction(-27))
        self.assertEqual(p.q, Fraction(search_integral.q_value(2, -3)))

    def test_rational_auxiliary_identity(self):
        p = search_rational.make_aux_point(1, -2, 2)

        self.assertEqual(p.q, p.y * p.y + p.x * p.y - p.x**3)
        self.assertEqual(p.x, Fraction(13, 2))
        self.assertEqual(p.y, Fraction(1, 2))

    def test_parse_z_values_deduplicates_and_rejects_nonpositive(self):
        self.assertEqual(search_rational.parse_z_values("1,2,2,3"), [1, 2, 3])

        with self.assertRaises(argparse.ArgumentTypeError):
            search_rational.parse_z_values("1,0")


if __name__ == "__main__":
    unittest.main()

## Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
## (bluescarni@gmail.com)
##
## This file is part of the pykep library.
##
## This Source Code Form is subject to the terms of the Mozilla
## Public License v. 2.0. If a copy of the MPL was not distributed
## with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

import unittest as _ut


def float_abs_error(a: float, b: float):
    return abs(a - b)

class anomaly_conversions_tests(_ut.TestCase):
    def test_m2e(self):
        import pykep as pk

        self.assertTrue(float_abs_error(pk.m2e(pk.e2m(0.1, 0.1), 0.1), 0.1) < 1e-14)

    def test_m2f(self):
        import pykep as pk

        self.assertTrue(float_abs_error(pk.m2f(pk.f2m(0.1, 0.1), 0.1), 0.1) < 1e-14)

    def test_f2e(self):
        import pykep as pk

        self.assertTrue(float_abs_error(pk.f2e(pk.e2f(0.1, 0.1), 0.1), 0.1) < 1e-14)

    def test_n2h(self):
        import pykep as pk

        self.assertTrue(float_abs_error(pk.n2h(pk.h2n(0.1, 10.1), 10.1), 0.1) < 1e-14)

    def test_n2f(self):
        import pykep as pk

        self.assertTrue(float_abs_error(pk.n2f(pk.f2n(0.1, 10.1), 10.1), 0.1) < 1e-14)

    def test_f2h(self):
        import pykep as pk

        self.assertTrue(float_abs_error(pk.f2h(pk.h2f(0.1, 10.1), 10.1), 0.1) < 1e-14)

    def test_f2zeta(self):
        import pykep as pk

        self.assertTrue(float_abs_error(pk.f2zeta(pk.zeta2f(0.1, 10.1), 10.1), 0.1) < 1e-14)


def run_test_suite():
    suite = _ut.TestSuite()
    suite.addTest(anomaly_conversions_tests("test_m2e"))
    suite.addTest(anomaly_conversions_tests("test_m2f"))
    suite.addTest(anomaly_conversions_tests("test_f2e"))
    suite.addTest(anomaly_conversions_tests("test_n2h"))
    suite.addTest(anomaly_conversions_tests("test_n2f"))
    suite.addTest(anomaly_conversions_tests("test_f2h"))
    suite.addTest(anomaly_conversions_tests("test_f2zeta"))

    test_result = _ut.TextTestRunner(verbosity=2).run(suite)

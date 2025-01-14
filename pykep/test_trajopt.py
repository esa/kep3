## Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
## (bluescarni@gmail.com)
##
## This file is part of the pykep library.
##
## This Source Code Form is subject to the terms of the Mozilla
## Public License v. 2.0. If a copy of the MPL was not distributed
## with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

import pykep as pk
import pygmo as pg
import numpy as np

import unittest as _ut


def float_rel_error(a: float, b: float):
    return abs(a - b) / abs(a)


class trajopt_mga_tests(_ut.TestCase):
    def test_construction(self):
        import pykep as pk

        earth = pk.planet(pk.udpla.jpl_lp("earth"))
        venus = pk.planet(pk.udpla.jpl_lp("venus"))
        udp = pk.trajopt.mga(
            seq=[earth, venus, earth, venus, earth],
            tof_encoding="direct",
            t0=[0, 1000],
            tof=[[30, 200], [30, 300], [30, 300], [30, 300]],
            vinf=2.5,
        )
        prob = pg.problem(udp)
        pop = pg.population(prob, 100)

    def test_encoding_to_encoding(self):
        import pykep as pk

        udp_direct = pk.trajopt.mga(tof_encoding="direct", tof=[[30, 200], [200, 300]])
        udp_alpha = pk.trajopt.mga(tof_encoding="alpha", tof=[230, 500])
        udp_eta = pk.trajopt.mga(tof_encoding="eta", tof=500)
        prob = pg.problem(udp_direct)
        pop = pg.population(prob, 100)
        x_direct = pop.champion_x
        gt = udp_direct.fitness(x_direct)[0]
        x_alpha = udp_alpha.direct2alpha(x_direct)
        x_eta = udp_eta.direct2eta(x_direct)
        self.assertTrue(float_rel_error(gt, udp_alpha.fitness(x_alpha)[0]) < 1e-14)
        self.assertTrue(float_rel_error(gt, udp_eta.fitness(x_eta)[0]) < 1e-14)
        self.assertTrue(
            float_rel_error(gt, udp_direct.fitness(udp_direct.alpha2direct(x_alpha))[0])
            < 1e-14
        )
        self.assertTrue(
            float_rel_error(gt, udp_direct.fitness(udp_eta.eta2direct(x_eta))[0])
            < 1e-14
        )


class gym_tests(_ut.TestCase):
    def cassini1(self):
        import pykep as pk

        udp = pk.trajopt.gym.cassini1
        # Ground truth checked by the old pykep code (up to 8 digits only as per differences with constants and all)
        x = [
            -554.5189290104555,
            103.27184879471751,
            335.41655259663474,
            80.50258543604521,
            862.0950563689543,
            2865.018040480413,
        ]
        f = udp.fitness(x)[0]
        self.assertTrue(float_rel_error(f, 80400.08898184073) < 1e-14)

        x = [
            -932.0532394941108,
            37.534681289972674,
            162.5093144821548,
            336.970139545233,
            1743.2915882004586,
            2527.8785180526465,
        ]
        f = udp.fitness(x)[0]
        self.assertTrue(float_rel_error(f, 217105.875031613573) < 1e-14)

        x = [
            -583.0776694390058,
            388.65047998036107,
            334.9959782156864,
            65.57508619540917,
            1520.2982946551908,
            2132.7771932619144,
        ]
        f = udp.fitness(x)[0]
        self.assertTrue(float_rel_error(f, 107218.08496509642) < 1e-14)

    def cassini2(self):
        import pykep as pk

        udp = pk.trajopt.gym.cassini2
        # Ground truth checked by the old pykep code (up to 8 digits only as per differences with constants and all)
        x = [-7.75699976e+02,  9.15777367e-01,  4.06442043e-01,  3.21309562e+03,
        6.81118341e-01,  1.62660490e+02, -1.58051063e+00,  1.28479507e+00,
        4.72699902e-01,  4.24319550e+02,  4.30475919e+00,  1.15739933e+00,
        2.55718252e-01,  5.44489098e+01, -1.54332794e+00,  1.27160729e+00,
        9.00000000e-01,  5.88481599e+02,  4.76774269e+00,  7.00000000e+01,
        1.00000000e-02,  2.20000000e+03]
        
        f = udp.fitness(x)[0]
        self.assertTrue(float_rel_error(f, 1511.7317645968126) < 1e-13)
        
    def rosetta(self):
        import pykep as pk

        udp = pk.trajopt.gym.rosetta
        # Ground truth checked by the old pykep code (up to 8 digits only as per differences with constants and all)
        x = [1.53488329e+03, 4.56388378e-01, 9.51717655e-01, 4.18212047e+03,
             4.32159299e-01, 3.65256539e+02, 5.03363275e+00, 2.38949977e+00,
             4.55746823e-01, 7.09999954e+02, 1.79894273e+00, 1.05000003e+00,
             6.09083347e-01, 2.60816142e+02, 4.95158968e+00, 3.16049580e+00,
             6.89049263e-01, 7.29775762e+02, 4.30823655e+00, 1.10842692e+00,
             4.16075410e-01, 1.84999995e+03]
        
        f = udp.fitness(x)[0]
        self.assertTrue(float_rel_error(f, 1371.4992633334382) < 1e-13)
        
    def eve_mga1dsm(self):
        import pykep as pk

        udp = pk.trajopt.gym.eve_mga1dsm
        # Ground truth checked by the old pykep code (up to 8 digits only as per differences with constants and all)
        x = [ 7.31864730e+02,  6.62420829e-01,  3.46714249e-01,  1.60589872e+03,
              4.69582854e-01,  2.98596210e+02, -1.90774304e+00,  2.05710242e+01,
              3.63127164e-01,  9.95555392e+01]
        f = udp.fitness(x)[0]
        self.assertTrue(float_rel_error(f, 47456.939061940415) < 1e-13)
        
    def eve_mga1dsm_a(self):
        import pykep as pk

        udp = pk.trajopt.gym.eve_mga1dsm_a
        # Ground truth checked by the old pykep code (up to 8 digits only as per differences with constants and all)
        x = [ 1.12861301e+02,  5.49233732e-01,  3.04597487e-02,  1.92554472e+03,
              5.22619135e-01,  8.46696560e-01, -2.64317289e+00,  2.16924824e+01,
              6.62441172e-01,  8.89339812e-02,  5.15383086e+02]
        f = udp.fitness(x)[0]
        self.assertTrue(float_rel_error(f, 1101622.7179572878) < 1e-13)
        
    def eve_mga1dsm_n(self):
        import pykep as pk

        udp = pk.trajopt.gym.eve_mga1dsm_n
        # Ground truth checked by the old pykep code (up to 8 digits only as per differences with constants and all)
        x = [5.86500918e+02, 6.32855532e-01, 9.59033298e-01, 1.93800759e+03,
             8.02901287e-01, 9.88679911e-01, 5.18276555e+00, 1.04655908e+01,
             5.84524787e-01, 9.68549775e-01]
        f = udp.fitness(x)[0]
        self.assertTrue(float_rel_error(f, 1917650.9004062244) < 1e-13)
        

class trajopt_mga1dsm_tests(_ut.TestCase):
    def test_construction(self):
        import pykep as pk

        earth = pk.planet(pk.udpla.jpl_lp("earth"))
        venus = pk.planet(pk.udpla.jpl_lp("venus"))
        udp = pk.trajopt.mga_1dsm(
            seq=[
                earth,
                venus,
                earth,
            ],
            t0=[0, 1000],
            tof=[[30, 200], [200, 300]],
            vinf=[0.5, 2.5],
            add_vinf_dep=False,
            add_vinf_arr=True,
            tof_encoding="direct",
            multi_objective=False,
            orbit_insertion=False,
            e_target=None,
            rp_target=None,
            eta_bounds=[0.1, 0.9],
            rp_ub=30,
        )
        prob = pg.problem(udp)
        pop = pg.population(prob, 100)

    def test_encoding_to_encoding(self):
        import pykep as pk

        udp_direct = pk.trajopt.mga_1dsm(
            tof_encoding="direct", tof=[[30, 200], [200, 300]]
        )
        udp_alpha = pk.trajopt.mga_1dsm(tof_encoding="alpha", tof=[230, 500])
        udp_eta = pk.trajopt.mga_1dsm(tof_encoding="eta", tof=500)
        prob = pg.problem(udp_direct)
        pop = pg.population(prob, 100)
        x_direct = pop.champion_x
        gt = udp_direct.fitness(x_direct)[0]
        x_alpha = udp_alpha.direct2alpha(x_direct)
        x_eta = udp_eta.direct2eta(x_direct, 500)
        self.assertTrue(float_rel_error(gt, udp_alpha.fitness(x_alpha)[0]) < 1e-14)
        self.assertTrue(float_rel_error(gt, udp_eta.fitness(x_eta)[0]) < 1e-14)
        self.assertTrue(
            float_rel_error(gt, udp_direct.fitness(udp_direct.alpha2direct(x_alpha))[0])
            < 1e-14
        )
        self.assertTrue(
            float_rel_error(gt, udp_direct.fitness(udp_eta.eta2direct(x_eta, 500))[0])
            < 1e-13
        )
        self.assertTrue(
            np.linalg.norm(udp_direct.alpha2direct(x_alpha) - x_direct) < 1e-13
        )
        self.assertTrue(
            np.linalg.norm(udp_direct.eta2direct(x_eta, 500) - x_direct) < 1e-13
        )

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


class gym_cassini1_tests(_ut.TestCase):
    def test_fitness(self):
        import pykep as pk

        udp = pk.trajopt.gym.cassini1
        # Three random values. Ground truth provided by the old pykep code
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
            < 1e-14
        )
        self.assertTrue(
            np.linalg.norm(udp_direct.alpha2direct(x_alpha) - x_direct) < 1e-14
        )
        self.assertTrue(
            np.linalg.norm(udp_direct.eta2direct(x_eta, 500) - x_direct) < 1e-14
        )

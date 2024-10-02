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

import unittest as _ut

def float_rel_error(a: float, b: float):
    return abs(a - b) / abs(a)

class trajopt_mga_tests(_ut.TestCase):
    def test_construction(self):
        import pykep as pk
        earth = pk.planet(pk.udpla.jpl_lp("earth"))
        venus = pk.planet(pk.udpla.jpl_lp("venus"))
        udp = pk.trajopt.mga(
            seq=[
                earth,
                venus,
                earth,
                venus,
                earth
            ],
            tof_encoding = "direct",
            t0=[0, 1000],
            tof=[[30, 200], [30, 300], [30, 300], [30, 300]],
            vinf=2.5,
        )
        prob = pg.problem(udp)
        pop = pg.population(prob, 100)

    def test_encoding_to_encoding(self):
        import pykep as pk
        udp_direct = pk.trajopt.mga(tof_encoding="direct", tof = [[30, 200], [200, 300]])
        udp_alpha = pk.trajopt.mga(tof_encoding="alpha", tof = [230, 500])
        udp_eta = pk.trajopt.mga(tof_encoding="eta", tof = 500)
        prob = pg.problem(udp_direct)
        pop = pg.population(prob, 100)
        x_direct = pop.champion_x
        gt = udp_direct.fitness(x_direct)[0]
        x_alpha = udp_alpha.direct2alpha(x_direct)
        x_eta = udp_eta.direct2eta(x_direct)
        self.assertTrue(float_rel_error(gt, udp_alpha.fitness(x_alpha)[0]) < 1e-14)
        self.assertTrue(float_rel_error(gt, udp_eta.fitness(x_eta)[0]) < 1e-14)
        self.assertTrue(float_rel_error(gt, udp_direct.fitness(udp_direct.alpha2direct(x_alpha))[0]) < 1e-14)
        self.assertTrue(float_rel_error(gt, udp_direct.fitness(udp_eta.eta2direct(x_eta))[0]) < 1e-14)

class gym_cassini1_tests(_ut.TestCase):
    def test_fitness(self):
        import pykep as pk
        udp = pk.trajopt.gym.cassini1
        # Three random values. Ground truth provided by a snapshot of the kep3 code validated agains the old pykep code 
        # (not identical, just validated ... as a the old pykep had a buggy jpl_lp eph method nhere fixed)
        x = [-554.5189290104555, 103.27184879471751, 335.41655259663474, 80.50258543604521, 862.0950563689543, 2865.018040480413  ]
        f = udp.fitness(x)[0]
        self.assertTrue(float_rel_error(f,80747.26667221037) < 1e-14)
        
        x = [-932.0532394941108, 37.534681289972674, 162.5093144821548, 336.970139545233, 1743.2915882004586, 2527.8785180526465  ]
        f = udp.fitness(x)[0]
        self.assertTrue(float_rel_error(f,216694.84791232392) < 1e-14)
        
        x = [-583.0776694390058, 388.65047998036107, 334.9959782156864, 65.57508619540917, 1520.2982946551908, 2132.7771932619144 ]
        f = udp.fitness(x)[0]
        self.assertTrue(float_rel_error(f,107787.63141728108) < 1e-14)
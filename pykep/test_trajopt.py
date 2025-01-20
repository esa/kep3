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
        self.assertTrue(float_rel_error(gt, udp_alpha.fitness(x_alpha)[0]) < 1e-13)
        self.assertTrue(float_rel_error(gt, udp_eta.fitness(x_eta)[0]) < 1e-13)
        self.assertTrue(
            float_rel_error(gt, udp_direct.fitness(udp_direct.alpha2direct(x_alpha))[0])
            < 1e-13
        )
        self.assertTrue(
            float_rel_error(gt, udp_direct.fitness(udp_eta.eta2direct(x_eta))[0])
            < 1e-13
        )


class gym_tests(_ut.TestCase):
    def test_cassini1(self):
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

    def test_cassini2(self):
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
        
    def test_rosetta(self):
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
        
    def test_eve_mga1dsm(self):
        import pykep as pk

        udp = pk.trajopt.gym.eve_mga1dsm
        # Ground truth checked by the old pykep code (up to 8 digits only as per differences with constants and all)
        x = [ 7.31864730e+02,  6.62420829e-01,  3.46714249e-01,  1.60589872e+03,
              4.69582854e-01,  2.98596210e+02, -1.90774304e+00,  2.05710242e+01,
              3.63127164e-01,  9.95555392e+01]
        f = udp.fitness(x)[0]
        self.assertTrue(float_rel_error(f, 47456.939061940415) < 1e-13)
        
    def test_eve_mga1dsm_a(self):
        import pykep as pk

        udp = pk.trajopt.gym.eve_mga1dsm_a
        # Ground truth checked by the old pykep code (up to 8 digits only as per differences with constants and all)
        x = [ 1.12861301e+02,  5.49233732e-01,  3.04597487e-02,  1.92554472e+03,
              5.22619135e-01,  8.46696560e-01, -2.64317289e+00,  2.16924824e+01,
              6.62441172e-01,  8.89339812e-02,  5.15383086e+02]
        f = udp.fitness(x)[0]
        self.assertTrue(float_rel_error(f, 1101622.7179572878) < 1e-13)
        
    def test_eve_mga1dsm_n(self):
        import pykep as pk

        udp = pk.trajopt.gym.eve_mga1dsm_n
        # Ground truth checked by the old pykep code (up to 8 digits only as per differences with constants and all)
        x = [5.86500918e+02, 6.32855532e-01, 9.59033298e-01, 1.93800759e+03,
             8.02901287e-01, 9.88679911e-01, 5.18276555e+00, 1.04655908e+01,
             5.84524787e-01, 9.68549775e-01]
        f = udp.fitness(x)[0]
        self.assertTrue(float_rel_error(f, 1917650.9004062244) < 1e-13)
        
    def test_juice(self):
        import pykep as pk

        udp = pk.trajopt.gym.juice
        # Ground truth checked by the old pykep code (up to 8 digits only as per differences with constants and all)
        x = [ 8.25587848945082e+03,  3.87927696743959e-01,
              5.27687480720340e-01,  2.06516101304215e+03,
              7.74309396597230e-01,  4.66381738991143e+02,
             -1.02751095385811e+00,  8.25250602985045e+00,
              1.85241667091363e-01,  1.10839235700056e+02,
              3.89240251780202e+00,  2.50359555867057e+00,
              6.46178864367641e-02,  2.15407037751815e+02,
              4.78686007014207e+00,  8.41382633987924e+00,
              3.38707892618604e-01,  3.19871173483077e+01,
              2.06974341215216e+00,  4.36373629930523e+00,
              4.66997711732296e-01,  7.22372043742636e+02,
              4.77835192401833e+00,  7.65391290702327e+00,
              4.89771110016942e-01,  1.03158288517734e+03]
        
        f = udp.fitness(x)[0]
        self.assertTrue(float_rel_error(f, 204.38446462495546) < 1e-13)
        
        # Solution previously in the old pykep gym code
        x = [ 8.16283083e+03,  6.41922787e-01,  6.51202691e-01,  2.51009414e+03,
        2.97841478e-01,  3.81541370e+02,  9.58572190e-01,  1.53007674e+00,
        3.06125365e-01,  1.49264351e+02,  4.10103546e+00,  2.39297670e+00,
        4.34424957e-01,  3.16066418e+02,  4.33225338e+00,  1.30946367e+00,
        4.52048883e-01,  1.63208108e+02,  5.93850330e-01,  1.34871269e+00,
        2.03288502e-01,  6.52494606e+02, -1.37902374e+00,  1.55482534e+00,
        1.96917559e-01,  1.08471126e+03]
        
        f = udp.fitness(x)[0]
        self.assertTrue(float_rel_error(f, -7.987614927397559) < 1e-13)
        
    def test_juice_mo(self):
        import pykep as pk

        udp = pk.trajopt.gym.juice_mo
        # Ground truth checked by the old pykep code (up to 8 digits only as per differences with constants and all)
        x = [ 8.08458777709711e+03,  3.98862191419023e-02,
              1.13145851845895e-01,  2.54109146808688e+03,
              7.47321603643624e-01,  2.27778122138468e-02,
             -4.41107070349122e+00,  6.27932399248486e+00,
              7.64033667254947e-01,  2.94097311318708e-01,
              2.72683135802239e+00,  3.48208267655823e+00,
              9.73171666093481e-02,  5.02201106930738e-01,
              3.76552737505492e+00,  5.97748097683895e+00,
              5.76935286152452e-01,  1.52899955890002e-01,
             -1.54719587734854e+00,  8.37088373080571e+00,
              7.83924298935545e-01,  5.32290337033144e-01,
             -5.62028545284311e+00,  3.37199051531738e+00,
              9.13407577042211e-01,  2.61007823689634e-01,
              2.27004720198764e+03]
        
        f = udp.fitness(x)
        self.assertTrue(float_rel_error(f[0], 427.31998557200325) < 1e-13)
        self.assertTrue(float_rel_error(f[1], 2270.0472019876433) < 1e-13)
        
    def test_messenger(self):
        import pykep as pk

        udp = pk.trajopt.gym.messenger
        # Ground truth checked by the old pykep code (up to 8 digits only as per differences with constants and all)
        x = [ 2.03241398e+03,  6.40762059e-01,  6.63357785e-01,  4.04989271e+03,
        6.63732323e-01,  4.50068524e+02, -3.86553343e+00,  3.52631372e+00,
        5.57888828e-01,  2.24619580e+02, -4.45910441e+00,  1.22736521e+00,
        7.08063036e-01,  2.17965497e+02, -2.47894274e+00,  1.43586128e+00,
        5.88391838e-01,  2.62423586e+02, -2.40594385e-02,  2.45470457e+00,
        7.25370468e-01,  3.58067954e+02,  1.47192632e+00,  1.05000000e+00,
        9.02984391e-01,  5.38436770e+02]
        
        f = udp.fitness(x)
        self.assertTrue(float_rel_error(f[0], 5855.81434335236) < 1e-13)

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
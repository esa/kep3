## Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
## (bluescarni@gmail.com)
##
## This file is part of the pykep library.
##
## This Source Code Form is subject to the terms of the Mozilla
## Public License v. 2.0. If a copy of the MPL was not distributed
## with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

import unittest as _ut
from .test_trajopt import *

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

class epoch_test(_ut.TestCase):
    def test_epoch_construction(self):
        import pykep as pk
        import datetime

        ep1 = pk.epoch(0., pk.epoch.julian_type.MJD2000)
        ep2 = pk.epoch(0., pk.epoch.julian_type.MJD)
        ep3 = pk.epoch(0., pk.epoch.julian_type.JD)
        self.assertTrue(ep1.mjd2000 == 0.0)
        self.assertTrue(ep2.mjd == 0.0)
        self.assertTrue(ep3.jd == 0.0)

        ep_py = datetime.datetime(2000, 1, 1, 0, 0, 0, 0)
        ep4 = pk.epoch(ep_py)
        self.assertTrue(ep4.mjd2000 == 0.0)
        self.assertRaises(TypeError, pk.epoch, datetime.timedelta(1.2))

        self.assertTrue(pk.epoch("2000-01") == ep1)
        self.assertTrue(pk.epoch("2000-01-01") == ep1)
        self.assertTrue(pk.epoch("2000-01-01T00") == ep1)
        self.assertTrue(pk.epoch("2000-01-01T00:00") == ep1)
        self.assertTrue(pk.epoch("2000-01-01T00:00:00") == ep1)


    def test_epoch_operators(self):
        import pykep as pk
        import datetime
        ep1 = pk.epoch(0.)
        ep2 = pk.epoch(1.)
        dt = datetime.timedelta(days=1)
        self.assertTrue(ep1 + dt == ep1 + 1.)
        self.assertTrue(ep2 > ep1)
        self.assertTrue(ep2 >= ep1)
        self.assertTrue(ep1 < ep2)
        self.assertTrue(ep1 <= ep2)
        self.assertTrue(ep1 == ep2 - 1)
        self.assertTrue(ep2 == ep1 + 1)

class my_udpla:
    def eph(self, ep):
        return [[1.,0.,0.],[0.,1.,0.]]
    def get_name(self):
        return "Boh"
    def get_mu_central_body(self):
        return 1.

class my_udpla_malformed1:
    def ephs(self, ep):
        return [[1.,0.,0.],[0.,1.,0.]]
    
class my_udpla_with_optionals:
    def eph(self, ep):
        # Some weird return value just to make it not constant
        return [[ep,0.,0.],[0.,1.,0.]]
    def eph_v(self, eps):
        retval = []
        [retval.extend([item,0.,0,.0,1.,.0]) for item in eps]
        return retval
    def elements(self, ep, el_type ):
        return [1.,2.,3.,4.,5.,6.]
    def period(self, mjd2000):
        return 3.14
    
class planet_test(_ut.TestCase):
    def test_planet_construction(self):
        import pykep as pk

        udpla_cpp = pk.udpla.keplerian(when = pk.epoch(0), elem = [1.,0.,0.,0.,0.,0.], mu_central_body = 1.)
        udpla_py = my_udpla()
        pla1 = pk.planet(udpla_cpp)
        pla2 = pk.planet(udpla_py)
        # Cannot construct from pk.planet
        self.assertRaises(TypeError, pk.planet, pla1)
        # Cannot construct from instance
        self.assertRaises(TypeError, pk.planet, pk.planet)
        # Cannot construct from malformed udpla
        self.assertRaises(NotImplementedError, pk.planet, my_udpla_malformed1())

    def test_udpla_extraction(self):
        import pykep as pk

        udpla_cpp = pk.udpla.keplerian(when = pk.epoch(0), elem = [1.,0.,0.,0.,0.,0.], mu_central_body = 1.)
        udpla_py = my_udpla()
        pla1 = pk.planet(udpla_cpp)
        pla2 = pk.planet(udpla_py)
        #self.assertTrue(pla1.extract(pk.udpla.keplerian) is udpla_cpp) questo fails ... perche'?
        self.assertTrue(pla2.extract(my_udpla) is udpla_py)
        self.assertTrue(pla1.extract(my_udpla) is None)
        self.assertTrue(pla2.extract(pk.udpla.keplerian) is None)
        self.assertTrue(pla1.is_(pk.udpla.keplerian))
        self.assertTrue(pla2.is_(my_udpla))
        self.assertTrue(not pla1.is_(my_udpla))
        self.assertTrue(not pla2.is_(pk.udpla.keplerian))

    def test_udpla_getters(self):
        import pykep as pk

        udpla_cpp = pk.udpla.keplerian(when = pk.epoch(0), elem = [1.,0.,0.,0.,0.,0.], mu_central_body = 1.56, added_params=[3., 2., 2.1])
        udpla_py = my_udpla()
        pla1 = pk.planet(udpla_cpp)
        pla2 = pk.planet(udpla_py)
        self.assertTrue(pla1.get_mu_central_body() == 1.56)
        self.assertTrue(pla1.get_mu_self() == 3.)
        self.assertTrue(pla1.get_radius() == 2.)
        self.assertTrue(pla1.get_safe_radius() == 2.1)
        self.assertTrue(pla2.get_mu_central_body() == 1.)
        self.assertTrue(pla2.get_mu_self() == -1)
        self.assertTrue(pla2.get_radius() == -1)
        self.assertTrue(pla2.get_safe_radius() == -1)

    def test_udpla_optional_methods(self):
        import pykep as pk
        import numpy as np
        # Testing eph_v
        udpla = my_udpla_with_optionals()
        pla = pk.planet(udpla)
        self.assertTrue(pla.elements(0.) == [1.,2.,3.,4.,5.,6.])
        r0, v0 = pla.eph(0.)
        r1, v1 = pla.eph(1.)
        self.assertTrue(np.all(pla.eph_v([0., 1]) == [r0+v0,r1+v1]))
        # Testing period
        self.assertTrue(pla.period() == 3.14)
        self.assertTrue(pla.period(when = 0.) == 3.14)
        self.assertTrue(pla.period(when = pk.epoch(0.)) == 3.14)
        # Testing elements
        self.assertTrue(pla.elements() == [1.,2.,3.,4.,5.,6.])
        self.assertTrue(pla.elements(when = 0.) == [1.,2.,3.,4.,5.,6.])
        self.assertTrue(pla.elements(when = pk.epoch(0.)) == [1.,2.,3.,4.,5.,6.])

class py_udplas_test(_ut.TestCase):
    def test_tle(self):
        import pykep as pk
        from sgp4.api import Satrec
        from sgp4 import exporter
        import numpy as np
        from pathlib import Path

        pk_path = Path(pk.__path__[0])
        data_file = str(pk_path / "data" / "tle.txt")

        ## Test eph
        file = open(data_file, "r")
        while(True):
            header = file.readline()
            if header == "":
                break
            line1 = file.readline()
            line2 = file.readline()
            udpla = pk.udpla.tle(line1 = line1, line2 = line2)
            pla = pk.planet(udpla)
            ref_epoch = pk.epoch("2023-10")
            rpk,vpk = pla.eph(ref_epoch)
            satellite = Satrec.twoline2rv(line1, line2)
            jd = ref_epoch.mjd2000 + 2451544.5
            jd_i = int(jd)
            jd_fr = jd-jd_i
            e, r, v = satellite.sgp4(jd_i, jd_fr)
            self.assertTrue(np.allclose(np.array(r) * 1000, rpk, equal_nan=True, atol=1e-13))
            self.assertTrue(np.allclose(np.array(v) * 1000, vpk, equal_nan=True, atol=1e-13))
        file.close()

        ## Test eph_v
        file = open(data_file, "r")
        while(True):
            header = file.readline()
            if header == "":
                break
            line1 = file.readline()
            line2 = file.readline()
            udpla = pk.udpla.tle(line1 = line1, line2 = line2)
            pla = pk.planet(udpla)
            ref_epoch = pk.epoch("2023-10")
            mjd2000s = np.linspace(ref_epoch.mjd2000, ref_epoch.mjd2000 + 10, 10)
            respk = pla.eph_v(mjd2000s)
            satellite = Satrec.twoline2rv(line1, line2)
            jds = [mjd2000 + 2451544.5 for mjd2000 in mjd2000s]
            jd_is = [int(item) for item in jds]
            jd_frs = [a-b for a,b in zip(jds, jd_is)]
            e, r, v = satellite.sgp4_array(np.array(jd_is), np.array(jd_frs))
            rv = np.hstack((r,v))
            rv = rv.reshape((-1,6))*1000
            self.assertTrue(np.allclose(rv, respk, equal_nan=True, atol=1e-13))

    def test_spice(self):
        import pykep as pk
        import spiceypy as pyspice
        import numpy as np
        from pathlib import Path


        pk_path = Path(pk.__path__[0])
        kernel_file = str(pk_path / "data" / "de440s.bsp")
        pk.utils.load_spice_kernels(kernel_file)

        # We test eph
        udpla = pk.udpla.spice("JUPITER BARYCENTER", "ECLIPJ2000", "SSB")
        pla = pk.planet(udpla)
        rvpk = udpla.eph(0.12345)
        rv, _ = pyspice.spkezr("JUPITER BARYCENTER", (0.12345-0.5)*pk.DAY2SEC, "ECLIPJ2000", "NONE", "SSB")
        self.assertTrue(np.allclose(rvpk[0], rv[0:3]*1000, atol=1e-13))
        self.assertTrue(np.allclose(rvpk[1], rv[3:]*1000, atol=1e-13))

        # We test eph_v
        mjd2000s = np.linspace(0.12345, 30, 100)
        rvpk = udpla.eph_v(mjd2000s)
        rv, _ = pyspice.spkezr("JUPITER BARYCENTER", (mjd2000s-0.5)*pk.DAY2SEC, "ECLIPJ2000", "NONE", "SSB")
        self.assertTrue(np.allclose(rvpk, np.array(rv)*1000, atol=1e-13))


class vsop2013_test(_ut.TestCase):
    def test_basic(self):
        import pykep as pk

        p = pk.planet(pk.udpla.vsop2013())
        self.assertTrue("mercury" in str(p))
        self.assertTrue("1e-05" in str(p))
        
class propagate_test(_ut.TestCase):
    def test_lagrangian(self):
        import pykep as pk
        import numpy as np
        
        r,v = pk.propagate_lagrangian(rv = [[1.23,0.12,-0.53],[0.0456,1.,0.2347623]], tof = 7.32, mu = 1.34, stm=False)
        r_gt = (0.04322525778936694, -1.3187148031966465, -0.3566980626026463)
        v_gt = (0.9223897138305034, 0.18875607625417432, -0.3722128151635408)
        self.assertTrue(np.allclose(r, r_gt, atol=1e-13))
        self.assertTrue(np.allclose(v, v_gt, atol=1e-13))

        
    def test_stark(self):
        import pykep as pk
        import numpy as np
        
        sp = pk.stark_problem(mu = 1.34)
        r_gt,v_gt = pk.propagate_lagrangian(rv = [[1.23,0.12,-0.53],[0.0456,1.,0.2347623]], tof = 7.32, mu = 1.34, stm=False)
        rvm = sp.propagate(rvm_state = [1.23,0.12,-0.53, 0.0456,1.,0.2347623, 1.], thrust = [0.,0.,0.], tof = 7.32)
        self.assertTrue(np.allclose(rvm[:3], r_gt, atol=1e-13))
        self.assertTrue(np.allclose(rvm[3:6], v_gt, atol=1e-13))
        
        rvm, _, _ = sp.propagate_var(rvm_state = [1.23,0.12,-0.53, 0.0456,1.,0.2347623, 100.], thrust = [1e-19,0.,0.], tof = 7.32)
        self.assertTrue(np.allclose(rvm[:3], r_gt, atol=1e-13))
        self.assertTrue(np.allclose(rvm[3:6], v_gt, atol=1e-13))

def compute_numerical_gradient(sf_leg, sf_leg_type = 'lf'):
    import numpy as np
    import pykep as pk
    import pygmo as pg

    state_length = np.array(sf_leg.rvs).flatten().size + 1
    throttle_length = np.array(sf_leg.throttles).size
    chromosome = np.zeros((state_length * 2 + throttle_length + 1))
    chromosome[0:state_length] = np.append(np.array(sf_leg.rvs).flatten(), sf_leg.ms)
    chromosome[state_length:state_length+throttle_length] = np.array(sf_leg.throttles)
    chromosome[state_length+throttle_length:state_length*2+throttle_length] = np.append(np.array(sf_leg.rvf).flatten(), sf_leg.mf)
    chromosome[-1] = sf_leg.tof

    def set_and_compute_constraints(chromosome, sf_leg_type = 'lf'):

        if sf_leg_type == 'hf' or sf_leg_type == 'high-fidelity':
            sf_leg_constraint = pk.leg.sims_flanagan_hf()
        else:
            sf_leg_constraint = pk.leg.sims_flanagan()
        sf_leg_constraint.cut = 0.5
        sf_leg_constraint.max_thrust = 1
        sf_leg_constraint.mu = 1
        sf_leg_constraint.isp = 1
        sf_leg_constraint.rvs = [chromosome[0:3],chromosome[3:6]]
        sf_leg_constraint.ms = chromosome[6]
        sf_leg_constraint.throttles = chromosome[state_length:state_length+throttle_length]
        sf_leg_constraint.rvf = [chromosome[state_length+throttle_length:state_length+throttle_length+3],chromosome[state_length+throttle_length+3:state_length+throttle_length+6]]
        sf_leg_constraint.mf = chromosome[2*state_length+throttle_length-1]
        sf_leg_constraint.tof = chromosome[2*state_length+throttle_length]
        eq_con = sf_leg_constraint.compute_mismatch_constraints()
        ineq_con = sf_leg_constraint.compute_throttle_constraints()
        return np.concatenate((eq_con, ineq_con))

    return pg.estimate_gradient_h(callable = lambda x : set_and_compute_constraints(x, sf_leg_type), x=chromosome)

class sims_flanagan_test(_ut.TestCase):


    def test_sims_flanagan(self):
        import numpy as np

        udpla_e = pk.udpla.vsop2013("earth_moon", 1e-2)
        udpla_j = pk.udpla.vsop2013("jupiter", 1e-2)
        earth = pk.planet(udpla_e)
        jupiter = pk.planet(udpla_j)
        dt_days = 1000
        dt = dt_days * pk.DAY2SEC
        t0 = 1233.3
        rv0 = earth.eph(t0)
        rv1 = jupiter.eph(t0 + dt_days)
        lp = pk.lambert_problem(rv0[0], rv1[0], dt, pk.MU_SUN)
        rv0[1] = lp.v0[0]
        rv1[1] = lp.v1[0]

        cut_values = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
        mc_list = []
        for i in range(1, 34):
            for cut in cut_values:
                throttles = [0.0] * i * 3
                sf_hf_leg = pk.leg.sims_flanagan(rv0, 1.0, throttles, rv1, 1.0, dt, 1.0, 1.0, pk.MU_SUN, cut)
                mc = sf_hf_leg.compute_mismatch_constraints()
                mc[0] /= pk.AU
                mc[1] /= pk.AU
                mc[2] /= pk.AU
                mc[3] /= pk.EARTH_VELOCITY
                mc[4] /= pk.EARTH_VELOCITY
                mc[5] /= pk.EARTH_VELOCITY
                mc[6] /= 1000
                mc_list.append(mc)
                
        self.assertTrue(np.array([np.max(i) < 1e-8 for i in mc_list]).all())

    def test_mc_grad(self):
        import numpy as np

        sf_leg = pk.leg.sims_flanagan()
        sf_leg.cut = 0.5
        sf_leg.throttles = np.array([0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24,
        0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34])
        sf_leg.rvs = np.array([[1, 0.1, -0.1], [0.2, 1.0, -0.2]])
        sf_leg.ms = 1
        sf_leg.rvf = np.array([[1.2, -0.1, 0.1], [-0.2, 1.023, -0.44]])
        sf_leg.mf = 13 / 15
        sf_leg.max_thrust = 1
        sf_leg.mu = 1
        sf_leg.isp = 1
        sf_leg.tof = 1
        #sf_leg.tol = 1e-16
        state_length = np.array(sf_leg.rvs).flatten().size + 1
        throttle_length = np.array(sf_leg.throttles).size

        num_grad = compute_numerical_gradient(sf_leg, sf_leg_type = 'lf')
        num_grad = num_grad.reshape((17, 45), order='C')
        grad_rvm, grad_rvm_bck, grad_final = sf_leg.compute_mc_grad()
        a_tc_grad = sf_leg.compute_tc_grad()
        a_grad = np.zeros((state_length+throttle_length // 3, 2 * state_length + throttle_length + 1))
        a_grad[0:state_length, 0:state_length] = grad_rvm
        a_grad[0:state_length, state_length:state_length + throttle_length] = grad_final[:,0:throttle_length] 
        a_grad[0:state_length, state_length+throttle_length:state_length*2+throttle_length] = grad_rvm_bck
        a_grad[0:state_length, state_length*2+throttle_length] = grad_final[:, throttle_length:throttle_length + 1].reshape(7,)
        a_grad[state_length:, state_length:state_length+throttle_length] = a_tc_grad
        self.assertTrue(np.allclose(num_grad, a_grad, atol=1e-8))

class sims_flanagan_hf_test(_ut.TestCase):
    def test_comparison_sf_and_sf_hf(self):
        import pykep as pk
        import numpy as np

        sf_leg = pk.leg.sims_flanagan()
        sf_leg.cut = 0.5
        sf_leg.throttles = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        sf_leg.rvs = np.array([[1, 0.1, -0.1], [0.2, 1.0, -0.2]])
        sf_leg.ms = 1
        sf_leg.rvf = np.array([[1.2, -0.1, 0.1], [-0.2, 1.023, -0.44]])
        sf_leg.mf = 13 / 15
        sf_leg.max_thrust = 1
        sf_leg.mu = 1
        sf_leg.isp = 1
        sf_leg.tof = 1
        rvm_mc_sf = sf_leg.compute_mismatch_constraints()

        sf_hf_leg = pk.leg.sims_flanagan_hf()
        sf_hf_leg.cut = 0.5
        sf_hf_leg.throttles = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        sf_hf_leg.rvms = np.array([1, 0.1, -0.1, 0.2, 1.0, -0.2, 1])
        sf_hf_leg.rvmf = np.array([1.2, -0.1, 0.1, -0.2, 1.023, -0.44, 13 / 15])
        sf_hf_leg.max_thrust = 1
        sf_hf_leg.mu = 1
        sf_hf_leg.isp = 1
        sf_hf_leg.tof = 1
        rvm_mc_sf_hf = sf_hf_leg.compute_mismatch_constraints()
        self.assertTrue(np.allclose(rvm_mc_sf, rvm_mc_sf_hf, atol=1e-13))

    def test_sims_flanagan_hf(self):
        import numpy as np

        udpla_e = pk.udpla.vsop2013("earth_moon", 1e-2)
        udpla_j = pk.udpla.vsop2013("jupiter", 1e-2)
        earth = pk.planet(udpla_e)
        jupiter = pk.planet(udpla_j)
        dt_days = 1000
        dt = dt_days * pk.DAY2SEC
        t0 = 1233.3
        rv0 = earth.eph(t0)
        rv1 = jupiter.eph(t0 + dt_days)
        lp = pk.lambert_problem(rv0[0], rv1[0], dt, pk.MU_SUN)
        rv0[1] = lp.v0[0]
        rv1[1] = lp.v1[0]


        cut_values = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
        mc_list = []
        for i in range(1, 34):
            for cut in cut_values:
                throttles = [0.0] * i * 3
                sf_hf_leg = pk.leg.sims_flanagan_hf(rv0, 1.0, throttles, rv1, 1.0, dt, 1.0, 1.0, pk.MU_SUN, cut)
                mc = sf_hf_leg.compute_mismatch_constraints()
                mc[0] /= pk.AU
                mc[1] /= pk.AU
                mc[2] /= pk.AU
                mc[3] /= pk.EARTH_VELOCITY
                mc[4] /= pk.EARTH_VELOCITY
                mc[5] /= pk.EARTH_VELOCITY
                mc[6] /= 1000
                mc_list.append(mc)
                
        self.assertTrue(np.array([np.max(i) < 1e-8 for i in mc_list]).all())
        
    def test_mc_grad_hf(self):
        import numpy as np

        sf_leg = pk.leg.sims_flanagan_hf()
        sf_leg.cut = 0.5
        sf_leg.throttles = np.array([0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24,
        0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34])
        sf_leg.rvs = np.array([[1, 0.1, -0.1], [0.2, 1.0, -0.2]])
        sf_leg.ms = 1
        sf_leg.rvf = np.array([[1.2, -0.1, 0.1], [-0.2, 1.023, -0.44]])
        sf_leg.mf = 13 / 15
        sf_leg.max_thrust = 1
        sf_leg.mu = 1
        sf_leg.isp = 1
        sf_leg.tof = 1
        state_length = np.array(sf_leg.rvs).flatten().size + 1
        throttle_length = np.array(sf_leg.throttles).size

        num_grad = compute_numerical_gradient(sf_leg, sf_leg_type = 'hf')
        num_grad = num_grad.reshape((17, 45), order='C')
        grad_rvm, grad_rvm_bck, grad_final = sf_leg.compute_mc_grad()
        a_tc_grad = sf_leg.compute_tc_grad()
        a_grad = np.zeros((state_length+throttle_length // 3, 2 * state_length + throttle_length + 1))
        a_grad[0:state_length, 0:state_length] = grad_rvm
        a_grad[0:state_length, state_length:state_length + throttle_length] = grad_final[:,0:throttle_length] 
        a_grad[0:state_length, state_length+throttle_length:state_length*2+throttle_length] = grad_rvm_bck
        a_grad[0:state_length, state_length*2+throttle_length] = grad_final[:, throttle_length:throttle_length + 1].reshape(7,)
        a_grad[state_length:, state_length:state_length+throttle_length] = a_tc_grad
        self.assertTrue(np.allclose(num_grad, a_grad, atol=1e-8))


def run_test_suite():
    tl = _ut.TestLoader()

    suite = _ut.TestSuite()
    suite.addTest(anomaly_conversions_tests("test_m2e"))
    suite.addTest(anomaly_conversions_tests("test_m2f"))
    suite.addTest(anomaly_conversions_tests("test_f2e"))
    suite.addTest(anomaly_conversions_tests("test_n2h"))
    suite.addTest(anomaly_conversions_tests("test_n2f"))
    suite.addTest(anomaly_conversions_tests("test_f2h"))
    suite.addTest(anomaly_conversions_tests("test_f2zeta"))
    suite.addTest(planet_test("test_planet_construction"))
    suite.addTest(planet_test("test_udpla_extraction"))
    suite.addTest(planet_test("test_udpla_getters"))
    suite.addTest(planet_test("test_udpla_optional_methods"))
    suite.addTest(epoch_test("test_epoch_construction"))
    suite.addTest(epoch_test("test_epoch_operators"))
    suite.addTest(propagate_test("test_lagrangian"))
    suite.addTest(propagate_test("test_stark"))
    suite.addTest(sims_flanagan_test("test_sims_flanagan"))
    suite.addTest(sims_flanagan_test("test_mc_grad"))
    suite.addTest(sims_flanagan_hf_test("test_comparison_sf_and_sf_hf"))
    suite.addTest(sims_flanagan_hf_test("test_sims_flanagan_hf"))
    suite.addTest(sims_flanagan_hf_test("test_mc_grad_hf"))
    suite.addTest(py_udplas_test("test_tle"))
    suite.addTest(py_udplas_test("test_spice"))
    suite.addTest(trajopt_mga_tests("test_construction"))
    suite.addTest(trajopt_mga_tests("test_encoding_to_encoding"))
    suite.addTest(gym_cassini1_tests("test_fitness"))

    suite.addTest(tl.loadTestsFromTestCase(vsop2013_test))



    test_result = _ut.TextTestRunner(verbosity=2).run(suite)

import pykep as pk
import numpy as np

import unittest as _ut


class zoh_test(_ut.TestCase):

    def test_zoh_mc(self):
        import numpy as np
        import pykep as pk

        # The integrators (Keplerian Propagation in Cartesian)
        tol=1e-12
        tol_var = 1e-6
        ta_global = pk.ta.get_zoh_kep(tol)
        ta_var_global = pk.ta.get_zoh_kep_var(tol_var)
        
        # We create a ballistic leg 
        t0 = 10000 # MJD2000
        t1 = 10400 # MJD2000
        pl0 = pk.planet(pk.udpla.jpl_lp("Earth"))
        pl1 = pk.planet(pk.udpla.jpl_lp("Mars"))
        r0, _ = pl0.eph(t0)
        r1, _ = pl1.eph(t1)
        l = pk.lambert_problem(r0=r0, r1=r1, tof = (t1-t0) * pk.DAY2SEC, mu = pk.MU_SUN)
        m0 = 1000
        m1 = 1000

        # nd units (in these units mu must be one as to use the tas)
        L = pk.AU
        MU = pk.MU_SUN
        TIME = np.sqrt(L**3/MU)
        V =  L/TIME
        ACC = V/TIME
        MASS = 1000
        F = MASS*ACC
        
        # We instantiate a ballistic leg and test for small mismatches
        n_trials = 50
        for i in range(n_trials):
            # leg data
            nseg = int(np.random.uniform(4, 20))
            veff = np.random.uniform(2000, 7000) * pk.G0

            # nd construction data
            state0 = [it/L for it in r0] + [it/V for it in l.v0[0]] + [m0/MASS]
            state1 = [it/L for it in r1] + [it/V for it in l.v1[0]] + [m1/MASS]
            veff_nd = veff / V
            tgrid = np.linspace(t0*pk.DAY2SEC/TIME, t1*pk.DAY2SEC/TIME, nseg+1)
            controls = np.random.uniform(-1,1, (4*nseg,)) * 1e-3
            controls[0::4] = 0 # zeroing the thrust magnitude
            cut = np.random.uniform(0,1)
            ta_global.pars[4] = 1. / veff_nd
            ta_var_global.pars[4] = 1. / veff_nd
            # construct the leg
            leg = pk.leg.zoh(state0, controls.tolist(), state1, tgrid, cut = cut, tas = [ta_global, ta_var_global])
            # test
            mc = leg.compute_mismatch_constraints()
            self.assertTrue(np.array([np.max(i) < 1e-12 for i in mc]).all())


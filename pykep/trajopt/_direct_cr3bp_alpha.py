## Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
## (bluescarni@gmail.com)## 
## This file is part of the kep3 library.## 
## This Source Code Form is subject to the terms of the Mozilla
## Public License v. 2.0. If a copy of the MPL was not distributed
## with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

import numpy as _np
import pykep as pk
import heyoka as _hy


class direct_cr3bp_alpha:
    """Represents the optimal low-thrust transfer between two cr3bp orbits (they can be propagated) using a direct method.

    This problem works internally using the :class:`~pykep.leg.sims_flanagan` and manipulates its initial and final states, as well as its transfer time T, final mass mf
    and the controls as to link the two orbits with a low-thrust trajectory.

    The decision vector is::

        z = [t0, mf, Vsx, Vsy, Vsz, Vfx, Vfy, Vfz, alphas, throttles, tof] - all in non-dimensional units (standard for CR3BP)

    where throttles is a vector of throttles structured as [u0x, u0y,u0z, ...]. By throttles we intend non dimensiona thrust levels in [0,1].
    and alphas is a vector of time-segments structured as [alpha0, alpha1, ...]. We can convert alphas to time segments using pk.alpha2direct(alphas, tof).

    """

    def __init__(
        self,
        pls=[[-0.1, 0.0, 0.0], [0.0, -0.1, 0.0]],
        plf=[[1, 0.0, 0.0], [0.0, +0.1, 0.0]],
        ms=1,
        mu=pk.MU_MOON/(pk.MU_EARTH+pk.MU_MOON),
        max_thrust=0.1,
        veff=30,
        t0_bounds=[0.0, 2*_np.pi],
        tof_bounds=[1.0, 2*_np.pi],
        mf_bounds=[1.0, 0.1],
        vinfs=0.1,
        vinff=0.0,
        nseg=10,
        cut=0.6,
        mass_scaling=1.0,
        r_scaling=1.0,
        v_scaling=1.0,
        with_gradient=True,
        high_fidelity=False,
        tas = (pk.ta.get_zero_hold_kep(1e-16),
                       pk.ta.get_zero_hold_kep_var(1e-16))
    ):
        """direct_cr3bp_alpha(pls, plf, ms = 1500, mu=_pk.MU_SUN, max_thrust=0.12, isp=3000, t0_bounds=[6700.0, 6800.0], tof_bounds=[200.0, 300.0], mf_bounds=[1300.0, 1500.0], vinfs=3.0, vinff=0.0, nseg=10, cut=0.6, mass_scaling=1500, r_scaling=pk.AU, v_scaling=pk.EARTH_VELOCITY, with_gradient=True)

        Args:
            *pls* (:class:`~pykep.planet`): Initial planet. Defaults to jpl_lp Earth.

            *plf* (:class:`~pykep.planet`): Final planet. Defaults to jpl_lp Mars.

            *ms* (:class:`float`): Initial spacecraft mass in kg. Defaults to 1000 kg.

            *mu* (:class:`float`): Gravitational parameter, default is for the Sun (:class:`~pykep.MU_SUN`).

            *max_thrust* (:class:`float`): Maximum thrust in Newtons. Defaults to 0.12 N.

            *veff* (:class:`float`): Exhaust velocity (Specific impulse * _pk.G0) in m/s. Defaults to 3000 s * _pk.G0.

            *t0_bounds* (:class:`list`): Bounds for departure epoch in MJD2000. Defaults to [6700.0, 6800.0].

            *tof_bounds* (:class:`list`): Bounds for time of flight in days. Defaults to [200, 300] days.

            *mf_bounds* (:class:`list`): Bounds for final mass in kg. Defaults to [1300.0, 1500.0] kg.

            *vinfs* (:class:`float`): Allowed magnitude for the departure's relative velocity in km/s. Defaults to 3.

            *vinff* (:class:`float`): Allowed magnitude for the arrival's relative velocity in km/s. Defaults to 0.

            *nseg* (:class:`int`): Number of segments for the trajectory. Defaults to 10.

            *cut* (:class:`float`): Cut parameter for the :class:`~pykep.leg.sims_flanagan`. Defaults to 0.6.

            *mass_scaling* (:class:`float`): Scaling factor for mass (used to scale constraints). Defaults to 1500.

            *r_scaling* (:class:`float`): Scaling factor for distance, (used to scale constraints). Defaults AU (:class:`~pykep.AU`).

            *v_scaling* (:class:`float`): Scaling factor for velocity (used to scale constraints). Defaults the Earth's velocity (:class:`~pykep.EARTH_VELOCITY`).

            *with_gradient* (:class:`bool`): Indicates if gradient information should be used. Defaults True.

            *high_fidelity* (:class:`bool`): Indicates if sims flanagan leg uses DV impulses, or zero-order hold continous thrust (note unclear how reliable graidents are). Defaults False.

        """
        # We add as data member one single Sims-Flanagan leg and set it using problem data
        if high_fidelity:
            self.leg = pk.leg.sims_flanagan_hf_alpha(
                tas = tas
            )
        else:
            raise("No sf_lf_alpha implemented for direct_cr3bp_alpha yet!")
        
        self.leg.ms = ms
        self.leg.max_thrust = max_thrust
        self.leg.veff = veff
        self.leg.mu = mu
        self.leg.cut = cut

        # We also need to propagate the target state
        self.ta_dyn = tas[0] # Extract from inputs
        self.ta_dyn.pars[0] = self.leg.mu
        self.ta_dyn.pars[1] = self.leg.veff # Veff
        self.ta_dyn.pars[2] = 0.0 # No thrust
        self.ta_dyn.pars[3] = 0.0 # No thrust
        self.ta_dyn.pars[4] = 0.0 # No thrust

        # Now make a cfunc for the dynamics
        x, y, z, vx, vy, vz, m = _hy.make_vars("x", "y", "z", "vx", "vy", "vz", "m")
        x_hy = _np.array([x, y, z, vx, vy, vz, m])
        self.cf_dyn = _hy.cfunc([row[1] for row in self.ta_dyn.sys], x_hy.tolist())

    
        # We define some additional datamembers useful later-on
        self.pls = pls
        self.plf = plf
        self.t0_bounds = t0_bounds
        self.tof_bounds = tof_bounds
        self.mf_bounds = mf_bounds
        self.vinfs = vinfs
        self.vinff = vinff
        self.nseg = nseg
        self.mass_scaling = mass_scaling
        self.r_scaling = r_scaling
        self.v_scaling = v_scaling
        if high_fidelity:
            self.thrust_scaling = 1 / max_thrust
        else:
            self.thrust_scaling = 1
        self.with_gradient = with_gradient
        self.high_fidelity = high_fidelity

    # z = [t0, mf, Vsx, Vsy, Vsz, Vfx, Vfy, Vfz, talphas, throttles, tof]
    def get_bounds(self):
        # Determine range on alphas
        factor = 0.25
        alpha_low = pk.direct2alpha([1/(self.nseg*factor)] * int(self.nseg*factor))[0][0]
        alpha_high = pk.direct2alpha([1/(self.nseg/factor)] * int(self.nseg/factor))[0][0]

        lb = (
            [self.t0_bounds[0], self.mf_bounds[0]]
            + [-self.vinfs] * 3  # in m/s.
            + [-self.vinff] * 3  # in m/s.
            + [alpha_low] * self.nseg
            + [-1, -1, -1] * self.nseg
            + [self.tof_bounds[0]]
        )
        ub = (
            [self.t0_bounds[1], self.mf_bounds[1]]
            + [self.vinfs] * 3  # in m/s.
            + [self.vinff] * 3  # in m/s.
            + [min(alpha_high,1-1e-3)] * self.nseg
            + [1, 1, 1] * self.nseg
            + [self.tof_bounds[1]]
        )
        return (lb, ub)

    def _set_leg_from_x(self, x):
        # We set the leg using data in the decision vector
        # Propagate initial states
        self.ta_dyn.time = 0.0
        self.ta_dyn.state[:3] = self.pls[0]
        self.ta_dyn.state[3:6] = self.pls[1]
        self.ta_dyn.propagate_until(x[0])
        rs = self.ta_dyn.state[:3].copy()
        vs = self.ta_dyn.state[3:6].copy()
        # Propagate target state
        self.ta_dyn.time = 0.0
        self.ta_dyn.state[:3] = self.plf[0]
        self.ta_dyn.state[3:6] = self.plf[1]
        self.ta_dyn.propagate_until(x[0] + x[-1])
        rf = self.ta_dyn.state[:3].copy()
        vf = self.ta_dyn.state[3:6].copy()
        # print('SF final state: ', rf, rs, x[-1])

        # We set the leg using data in the decision vector
        self.leg.rvs = [rs, [a + b for a, b in zip(vs, x[2:5])]]  # we add vinfs
        self.leg.rvf = [rf, [a + b for a, b in zip(vf, x[5:8])]]  # we add vinff
        self.leg.tof = x[-1] # * _pk.DAY2SEC
        self.leg.mf = x[1]

        # Split alphas and throttles
        data = x[8:-1]
        alphas = data[:self.nseg]
        throttles = data[self.nseg:]

        # Decode alphas to direct
        T = pk.alpha2direct(alphas, x[-1])
        # print('Alphas: ', alphas, 'T', T)

        # Now save the modified back to self.leg.talphas and self.leg.throttles
        self.leg.talphas = T
        self.leg.throttles = throttles

        # We return the eph as to avoid having to recompute them later on if needed.
        return rs, vs, rf, vf

    # z = [t0, mf, Vsx, Vsy, Vsz, Vfx, Vfy, Vfz, talphas, throttles, tof]
    def fitness(self, x):
        # 1 - We set the optimality principle
        mf = x[1]

        # 2 - We compute the constraints violations (mismatch+throttle)
        self._set_leg_from_x(x)  # set the leg
        ceq = self.leg.compute_mismatch_constraints()
        cineq = self.leg.compute_throttle_constraints()

        # 3 - We add the departure vinfs constraint (quadratic)
        cineq = cineq + [
            (x[2] ** 2 + x[3] ** 2 + x[4] ** 2 - self.vinfs**2) / (self.v_scaling**2)
        ] 
        # We add the departure vinff constraint (quadratic)
        cineq = cineq + [
            (x[5] ** 2 + x[6] ** 2 + x[7] ** 2 - self.vinff**2) / (self.v_scaling**2)
        ]

        # 4 - We assemble the return fitness
        retval = _np.array([-mf] + ceq + cineq)  # here we can sum lists

        # 5 - We scale the values in nd units (numerical solvers are sensitive to well-scaled values)
        retval[0] /= self.mass_scaling
        retval[1:4] /= self.r_scaling
        retval[4:7] /= self.v_scaling
        retval[7] /= self.mass_scaling

        return retval

    def has_gradient(self):
        return False

    def get_nec(self):
        return 7

    def get_nic(self):
        return self.nseg + 2

    def pretty(self, x):
        """
        Prints a human readable representation of the transfer.

        Args:
            *x* (:class:`list`): The decision vector containing final mass, thrust direction, and time of flight.
        """
        self._set_leg_from_x(x)
        print(f"\nLow-thrust NEP transfer")
        print(f"Departure: {self.pls}\nArrival: {self.plf}")
        print(
            f"\nLaunch epoch: {x[0]:.5f}"
        )
        print(
            f"Arrival epoch: {x[0]+x[-1]:.5f}"
        )

        print(f"Time of flight (-): {x[-1]:.5f} ")
        print(
            f"\nLaunch DV (-) {_np.sqrt(x[2] ** 2 + x[3] ** 2 + x[4] ** 2):.8f} - [{x[2]},{x[3]},{x[4]}]"
        )
        print(
            f"Arrival DV (-) {_np.sqrt(x[5] ** 2 + x[6] ** 2 + x[7] ** 2):.8f} - [{x[5]},{x[6]},{x[7]}]"
        )
        print(f"Final mass (mf/m0): {x[1]}")
        print(f"\nDetails on the low-thrust leg: (NEEDS fixing) ")
        print(self.leg)

    def plot(
        self,
        x,
        ax=None,
        units=1,
        show_midpoints=False,
        show_gridpoints=False,
        show_throttles=False,
        length=0.1,
        arrow_length_ratio=0.05,
        **kwargs,
    ):
        """
        Plots the trajectory leg  3D axes.

        Args:
            *x* (:class:`list`): The decision vector containing final mass, thrust direction, and time of flight.

            *ax* (:class:`mpl_toolkits.mplot3d.axes3d.Axes3D`, optional): The 3D axis to plot on. Defaults to None.

            *units* (:class:`float`, optional): The unit scale for the plot. Defaults to _pk.AU.

            *show_midpoints* (:class:`bool`, optional): Whether to show midpoints on the trajectory. Defaults to False.

            *show_gridpoints* (:class:`bool`, optional): Whether to show grid points on the trajectory. Defaults to False.

            *show_throttles* (:class:`bool`, optional): Whether to show throttle vectors. Defaults to False.

            *length* (:class:`float`, optional): Length of the throttle vectors. Defaults to 0.1.

            *arrow_length_ratio* (:class:`float`, optional): Arrow length ratio for the throttle vectors. Defaults to 0.05.

            *\\*\\*kwargs*: Additional keyword arguments for the plot.

        Returns:
            :class:`mpl_toolkits.mplot3d.axes3d.Axes3D`: The 3D axis with the plotted trajectory.
        """
        self._set_leg_from_x(x)
        sf = self.leg
        # Making the axis
        if ax is None:
            ax = pk.plot.make_3Daxis(figsize=(7, 7))
            
        rs, _ = sf.rvs
        rf, _ = sf.rvf

        # Propagate initial states
        tgrid = _np.linspace(0.0, 2*(x[0] + x[-1]), 1000)
        self.ta_dyn.time = 0.0
        self.ta_dyn.state[:3] = self.pls[0]
        self.ta_dyn.state[3:6] = self.pls[1]
        sol_0 = self.ta_dyn.propagate_grid(tgrid)
        status = sol_0[0]
        integration_0 = sol_0[5]
        # Crop time
        tgrid = tgrid[:len(integration_0[:,0])]
        if status.name != 'time_limit':
            integration_0 = _np.concatenate([integration_0, [self.ta_dyn.state]]) 

        # Propagate target state
        self.ta_dyn.time = 0.0
        self.ta_dyn.state[:3] = self.plf[0]
        self.ta_dyn.state[3:6] = self.plf[1]
        sol_f = self.ta_dyn.propagate_grid(tgrid)
        status = sol_f[0]
        integration_f = sol_f[5]
        # Crop time
        tgrid = tgrid[:len(integration_f[:,0])]
        if status.name != 'time_limit':
            integration_f = _np.concatenate([integration_f, [self.ta_dyn.state]]) 

        # Extract traj_x0ectory from integration results
        traj_x0 = integration_0[:, :3]  # x, y, z
        traj_xf = integration_f[:, :3]  # x, y, z

        # MATLAB tab blue color
        matlab_colors = [
            (0, 0.4470, 0.7410),
            (0.8500, 0.3250, 0.0980),
            (0.9290, 0.6940, 0.1250),
            (0.4940, 0.1840, 0.5560),
            (0.4660, 0.6740, 0.1880),
            (0.3010, 0.7450, 0.9330),
            (0.6350, 0.0780, 0.1840)
        ]

        ax.plot(traj_x0[:, 0], traj_x0[:, 1], traj_x0[:, 2], color=matlab_colors[0], label='x0')
        ax.plot(traj_xf[:, 0], traj_xf[:, 1], traj_xf[:, 2], color=matlab_colors[2], label='xf')


        # Plotting the trajctory leg
        if self.high_fidelity:
            ax = pk.plot.add_sf_hf_leg(
                ax,
                sf,
                units=units,
                show_throttles=show_throttles,
                length=length,
                show_gridpoints=show_gridpoints,
                arrow_length_ratio=arrow_length_ratio,
                use_alpha = True,
                N = 20,
                **kwargs,
            )
        else:
            raise("No sf_lf_alpha implemented for direct_cr3bp_alpha yet!")  
        return ax
import numpy as _np
import pykep as _pk
from copy import deepcopy as _deepcopy

_posvel0 = [
    [34110913367.783306, -139910016918.87585, -14037825669.025244],
    [29090.9902134693, 10000.390168313803, 1003.3858682643288],
]
_posvelf = [
    [-159018773159.22266, -18832495968.945133, 15781467087.350443],
    [2781.182556622003, -28898.40730995848, -483.4533989771214],
]


class pontryagin_cartesian_mass:
    """
    This class is a pygmo (http://esa.github.io/pygmo/) UDP representing a single-shooting indirect method
    for the minimum mass, fixed time, optimization of a low-thrust trajectory.

    The decision vector is::

      x = [lx, ly, lz, lvx, lvy, lvz, lm, l0]

    .. note::

       Gradients are also provided via the variational version of the augmented dynamics.
    """

    def __init__(
        self,
        posvel0=_posvel0,
        posvelf=_posvelf,
        tof=250,
        mu=_pk.MU_SUN,
        eps=1e-4,
        T_max=0.6,
        Isp=3000,
        m0=1500,
        L=_pk.AU,
        TIME=_pk.YEAR2DAY * _pk.DAY2SEC,
        MASS=1500,
        with_gradient=False,
    ):
        r"""pykep.trajopt.pontryagin_cartesian(start=default, target=default, tof=250, mu=1.32712440018e+20, eps=1e-4, T_max=0.6, Isp=3000, m0=1500, L=1.495978707e+11, TIME=31557600.0, MASS=1500, with_gradient=False)

        Args:
            *posvel0* (:class:`list`): the initial position and velocity of the spacecraft.

            *posvelf* (:class:`list`): the final position and velocity of the spacecraft.

            *tof* (:class:`float`): the time of flight of the trajectory (in days).

            *mu* (:class:`float`): the gravitational parameter of the central body.

            *eps* (:class:`float`): the accuracy of the numerical integrator.

            *T_max* (:class:`float`): the maximum thrust of the spacecraft.

            *Isp* (:class:`float`): the specific impulse of the spacecraft.

            *m0* (:class:`float`): the initial mass of the spacecraft.

            *L* (:class:`float`): units for length. All inputs will be scaled to these units and must be of the same order as to obtain a well scaled problem.

            *TIME* (:class:`float`): units for time. All inputs will be scaled to these units and should be of the same order as to obtain a well scaled problem.

            *MASS* (:class:`float`): units for mass. All inputs will be scaled to these units and should be of the same order as to obtain a well scaled problem.

            *with_gradient* (:class:`bool`): whether to use the gradient of the constraints or not.
        """
        # Non-dimensional units
        VEL = L / TIME  # Unit for velocity (1 AU/year)
        ACC = VEL / TIME  # Unit for acceleration (1 AU/year^2)

        # We redefine the user inputs in non dimensional units
        self.mu = mu / (L**3 / TIME**2)
        self.eps = eps
        self.c1 = T_max / (MASS * ACC)
        self.c2 = (T_max / Isp / _pk.G0) / MASS * TIME

        self.posvel0 = [[it / L for it in posvel0[0]], [it / VEL for it in posvel0[1]]]
        self.posvelf = [[it / L for it in posvelf[0]], [it / VEL for it in posvelf[1]]]

        self.m0 = m0 / MASS
        self.tof = tof * _pk.DAY2SEC / TIME

        self.ta = _deepcopy(_pk.ta.get_pc(1e-16, _pk.optimality_type.MASS))
        self.ta_var = _deepcopy(_pk.ta.get_pc_var(1e-8, _pk.optimality_type.MASS))
        self.ic_var = _deepcopy(self.ta_var.state[14:])

        self.MASS = MASS
        self.L = L
        self.TIME = TIME

        self.with_gradient = with_gradient

    def get_bounds(self):
        lb = [-1.0] * 7 + [0.0]
        ub = [1.0] * 8
        return [lb, ub]

    def set_ta_state(self, x):
        # Preparing the numerical integration parameters
        self.ta.pars[0] = self.mu
        self.ta.pars[1] = self.c1
        self.ta.pars[2] = self.c2
        self.ta.pars[3] = self.eps
        self.ta.pars[4] = x[7]
        self.ta.time = 0.0

        # And initial conditions
        self.ta.state[:3] = self.posvel0[0]
        self.ta.state[3:6] = self.posvel0[1]
        self.ta.state[6] = self.m0
        self.ta.state[7:14] = x[:7]

    def set_ta_var_state(self, x):
        # Preparing the numerical integration parameters
        self.ta_var.pars[0] = self.mu
        self.ta_var.pars[1] = self.c1
        self.ta_var.pars[2] = self.c2
        self.ta_var.pars[3] = self.eps
        self.ta_var.pars[4] = x[7]
        self.ta_var.time = 0.0

        # And initial conditions
        self.ta_var.state[:3] = self.posvel0[0]
        self.ta_var.state[3:6] = self.posvel0[1]
        self.ta_var.state[6] = self.m0
        self.ta_var.state[7:14] = x[:7]
        self.ta_var.state[14:] = self.ic_var

    def fitness(self, x):
        # Single Shooting
        self.set_ta_state(x)
        self.ta.propagate_until(self.tof)
        # Assembling the constraints
        ceq = [(self.ta.state[0] - self.posvelf[0][0])]
        ceq += [(self.ta.state[1] - self.posvelf[0][1])]
        ceq += [(self.ta.state[2] - self.posvelf[0][2])]
        ceq += [(self.ta.state[3] - self.posvelf[1][0])]
        ceq += [(self.ta.state[4] - self.posvelf[1][1])]
        ceq += [(self.ta.state[5] - self.posvelf[1][2])]
        ceq += [self.ta.state[13]]  # lm = 0
        ceq += [sum([it * it for it in x]) - 1.0]  # |lambdas|^2 = 1
        return [1.0] + ceq

    def gradient(self, x):
        # Single Shooting of variational equations
        self.set_ta_var_state(x)
        self.ta_var.propagate_until(self.tof)
        # Assembling the gradient
        # Fitness (sparsity takes care of the fitness since it does not depend on the dv)
        retval = []
        # Constraints
        for i in range(6):
            sl = self.ta_var.get_vslice(order=1, component=i)
            retval.extend(list(self.ta_var.state[sl]))
        # Constraint on mass
        sl = self.ta_var.get_vslice(order=1, component=13)
        retval.extend(list(self.ta_var.state[sl]))
        # Norm constraint
        for item in x:
            retval.extend([2 * item])
        return retval

    def gradient_sparsity(self):
        retval = []
        for i in range(1, 9):
            for j in range(8):
                retval.append((i, j))
        return retval

    def has_gradient(self):
        return self.with_gradient

    def get_nec(self):
        return 8

    def plot_trajectory(self, x, N=100, ax3D=None):
        """
        This function plots the trajectory encoded in the decision vector x.

        Args:
            *x* (:class:`list`): the decision vector.
            *N* (:class:`int`): the number of points to use in the plot.
            *ax3D* (:class:`matplotlib.axes._axes.Axes`): the axis to use for the plot. If None, a new axis is created.

        Returns:
            *ax3D* (:class:`matplotlib.axes._axes.Axes`): the axis of the plot.
        """
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d.art3d import Line3DCollection
        import matplotlib.colors as mcolors

        # Single Shooting
        self.set_ta_state(x)
        t_grid = _np.linspace(0, self.tof, N)
        sol = self.ta.propagate_grid(t_grid)
        # We make the axis if needed
        if ax3D is None:
            ax3D = _pk.plot.make_3Daxis()
        # Adding the main body
        _pk.plot.add_sun(ax3D)
        # Adding the osculating orbits to the initial conditions
        pl1 = _pk.planet(
            _pk.udpla.keplerian(
                when=_pk.epoch(0), posvel=self.posvel0, mu_central_body=self.mu
            )
        )
        pl2 = _pk.planet(
            _pk.udpla.keplerian(
                when=_pk.epoch(0), posvel=self.posvelf, mu_central_body=self.mu
            )
        )
        _pk.plot.add_planet_orbit(ax3D, pl1, c="gray", units=1)
        _pk.plot.add_planet_orbit(ax3D, pl2, c="gray", units=1)
        # Plotting the trajectory
        # Assuming sol[-1] is an array of shape (N, 3)
        state = sol[-1]  # Extract (x, y, z) coordinates
        x, y, z = state[:, 0], state[:, 1], state[:, 2]

        # Compute color values
        u_func = _pk.ta.get_pc_u_cfunc(_pk.optimality_type.MASS)
        c = u_func(
            _np.ascontiguousarray(state.T),
            pars=_np.ascontiguousarray(_np.tile(self.ta.pars, (N, 1)).T),
        )
        c = _np.ravel(c)  # Flatten if needed
        # Create segments for the curve
        points = _np.array([x, y, z]).T.reshape(-1, 1, 3)  # Shape: (N, 1, 3)
        segments = _np.concatenate(
            [points[:-1], points[1:]], axis=1
        )  # Shape: (N-1, 2, 3)

        # Normalize c values for colormap
        cmap = plt.get_cmap("coolwarm")  # Blue (low) to red (high)
        norm = mcolors.Normalize(vmin=_np.min(c), vmax=_np.max(c))

        # Create a Line3DCollection
        lc = Line3DCollection(segments, cmap=cmap, norm=norm, linewidth=2)
        lc.set_array(c[:-1])  # Assign colors to segments

        # Add the colored line to the plot
        ax3D.add_collection(lc)
        # Plotting the boundary conditions
        ax3D.scatter(self.posvel0[0][0], self.posvel0[0][1], self.posvel0[0][2])
        ax3D.scatter(self.posvelf[0][0], self.posvelf[0][1], self.posvelf[0][2])
        return ax3D

    def plot_misc(self, x, N=100):
        """
        This function plots the throttle, thrust direction, switching function, mass costate and Hamiltonian
        of the trajectory encoded in the decision vector x.

        Args:
            *x* (:class:`list`): the decision vector.
            *N* (:class:`int`): the number of points to use in the plot.

        Returns:
            *axs* (:class:`list`): the list of axis of the plots.
        """
        from matplotlib import pyplot as plt

        # Single Shooting
        self.set_ta_state(x)
        t_grid = _np.linspace(0, self.tof, N)
        sol = self.ta.propagate_grid(t_grid)
        # Retreive useful cfuncs
        # The Hamiltonian
        H_func = _pk.ta.get_pc_H_cfunc(_pk.optimality_type.MASS)
        # The switching function
        SF_func = _pk.ta.get_pc_SF_cfunc(_pk.optimality_type.MASS)
        # The magnitude of the throttle
        u_func = _pk.ta.get_pc_u_cfunc(_pk.optimality_type.MASS)
        # The thrust direction
        i_vers_func = _pk.ta.get_pc_i_vers_cfunc(_pk.optimality_type.MASS)
        # Create axis
        _, axs = plt.subplots(3, 2, figsize=(10, 10))
        axs[1, 0].set_title("Mass")
        axs[1, 0].plot(t_grid, self.MASS * sol[-1][:, 6].T)
        # Plot throttle
        throttle = u_func(
            _np.ascontiguousarray(sol[-1].T),
            pars=_np.ascontiguousarray(_np.tile(self.ta.pars, (N, 1)).T),
        )
        axs[0, 0].set_title("Throttle (log-scale)")
        axs[0, 0].semilogy(t_grid, _np.squeeze(throttle))
        # Plot thrust direction
        thrust_dir = i_vers_func(_np.ascontiguousarray(sol[-1][:, 10:13].T))
        for i in range(3):
            axs[0, 1].plot(t_grid, thrust_dir[i, :])
        axs[0, 1].set_title("Thrust direction")
        # Plot switching function
        SF = SF_func(
            _np.ascontiguousarray(sol[-1].T),
            pars=_np.ascontiguousarray(_np.tile(self.ta.pars, (N, 1)).T),
        )
        axs[1, 1].set_title("Switching Function")
        axs[1, 1].plot(t_grid, _np.squeeze(SF))
        axs[1, 1].hlines(0, t_grid[0], t_grid[-1], "k")
        # Plot mass costate
        axs[2, 0].set_title("lm")
        axs[2, 0].plot(sol[-1][:, 13].T)
        # Plot Hamiltonian
        Ham = H_func(
            _np.ascontiguousarray(sol[-1].T),
            pars=_np.ascontiguousarray(_np.tile(self.ta.pars, (N, 1)).T),
        )
        axs[2, 1].set_title("Hamiltonian")
        axs[2, 1].plot(t_grid, _np.squeeze(Ham))

        plt.tight_layout()

        return axs


# Adding the osculating orbits to the initial conditions
_pl0 = _pk.planet(
    _pk.udpla.keplerian(when=_pk.epoch(0), posvel=_posvel0, mu_central_body=_pk.MU_SUN)
)
_plf = _pk.planet(
    _pk.udpla.keplerian(
        when=_pk.epoch(0) + 250, posvel=_posvelf, mu_central_body=_pk.MU_SUN
    )
)


class pontryagin_cartesian_time:
    """
    This class is a pygmo (http://esa.github.io/pygmo/) UDP representing a single-shooting indirect method
    for the minimum time, optimization of a low-thrust trajectory.

    The decision vector is::

      x = [lx, ly, lz, lvx, lvy, lvz, lm, lJ, tof]

    .. note::

       Gradients are also provided via the variational version of the augmented dynamics.
    """

    def __init__(
        self,
        source=_pl0,
        target=_plf,
        t0=_pk.epoch(0),
        tof_guess=250,
        T_max=0.6,
        Isp=3000,
        m0=1500,
        L=_pk.AU,
        TIME=_pk.YEAR2DAY * _pk.DAY2SEC,
        MASS=1500,
        with_gradient=False,
    ):
        r"""pykep.trajopt.pontryagin_cartesian(start=default, target=default, tof=250, mu=1.32712440018e+20, eps=1e-4, T_max=0.6, Isp=3000, m0=1500, L=1.495978707e+11, TIME=31557600.0, MASS=1500, with_gradient=False)

        Args:
            *source* (:class:`list`): the initial planet.

            *target* (:class:`list`): the final planet.

            *t0* (:class:`~pk.epoch`): the initial epoch.

            *tof_guess* (:class:`float`): a guess for the time of flight. Bound will be defined as tof_guess/2 and 2*tof_guess.

            *T_max* (:class:`float`): the maximum thrust of the spacecraft.

            *Isp* (:class:`float`): the specific impulse of the spacecraft.

            *m0* (:class:`float`): the initial mass of the spacecraft.

            *L* (:class:`float`): units for length. All inputs will be scaled to these units and must be of the same order as to obtain a well scaled problem.

            *TIME* (:class:`float`): units for time. All inputs will be scaled to these units and should be of the same order as to obtain a well scaled problem.

            *MASS* (:class:`float`): units for mass. All inputs will be scaled to these units and should be of the same order as to obtain a well scaled problem.

            *with_gradient* (:class:`bool`): whether to use the gradient of the constraints or not.
        """
        # Non-dimensional units
        VEL = L / TIME  # Unit for velocity (1 AU/year)
        ACC = VEL / TIME  # Unit for acceleration (1 AU/year^2)
        
        # Store user inputs
        self.Isp = Isp
        self.t0 = t0

        # We redefine the user inputs in non dimensional units.
        self.mu = source.get_mu_central_body() / (L**3 / TIME**2)
        self.c1 = T_max / (MASS * ACC)
        self.c2 = (T_max / Isp / _pk.G0) / MASS * TIME
        self.m0 = m0 / MASS
        self.tof_guess = tof_guess * _pk.DAY2SEC / TIME

        # Initial position is computed once only upon construction.
        r0, v0 = source.eph(t0)
        r0 = [it / L for it in r0]
        v0 = [it / VEL for it in v0]
        self.posvel0 = [r0, v0]
        
        # Storing the planets
        self.source = source
        self.target = target

        # And the Taylor integrators
        self.ta = _deepcopy(_pk.ta.get_pc(1e-16, _pk.optimality_type.TIME))
        self.ta_var = _deepcopy(_pk.ta.get_pc_var(1e-8, _pk.optimality_type.TIME))
        self.ic_var = _deepcopy(self.ta_var.state[14:])

        # Compiled functions
        self.dyn_func = _pk.ta.get_pc_dyn_cfunc(_pk.optimality_type.TIME)

        # Non dimensional units
        self.MASS = MASS
        self.L = L
        self.TIME = TIME
        self.VEL = VEL
        self.ACC = ACC

        # Boolean
        self.with_gradient = with_gradient

    def get_bounds(self):
        lb = [-1.0] * 6 + [0.0] * 2 + [self.tof_guess / 2]
        ub = [1.0] * 8 + [self.tof_guess * 2]
        return [lb, ub]

    def set_ta_state(self, x):
        # Preparing the numerical integration parameters
        self.ta.pars[0] = self.mu
        self.ta.pars[1] = self.c1
        self.ta.pars[2] = self.c2
        self.ta.time = 0.0

        # And initial conditions
        self.ta.state[:3] = self.posvel0[0]
        self.ta.state[3:6] = self.posvel0[1]
        self.ta.state[6] = self.m0
        self.ta.state[7:14] = x[:7]

    def set_ta_var_state(self, x):
        # Preparing the variational numerical integration parameters
        self.ta_var.pars[0] = self.mu
        self.ta_var.pars[1] = self.c1
        self.ta_var.pars[2] = self.c2
        self.ta_var.time = 0.0

        # And initial conditions
        self.ta_var.state[:3] = self.posvel0[0]
        self.ta_var.state[3:6] = self.posvel0[1]
        self.ta_var.state[6] = self.m0
        self.ta_var.state[7:14] = x[:7]
        self.ta_var.state[14:] = self.ic_var

    def fitness(self, x):
        # x = [lx, ly, lz, lvx, lvy, lvz, lm, lJ, tof]
        # Single Shooting
        self.set_ta_state(x[:8])
        self.ta.propagate_until(x[8])
        
        # Computing the target position at epoch
        rf, vf = self.target.eph(self.t0 + x[8] * self.TIME * _pk.SEC2DAY)
        rf = [it / self.L for it in rf]
        vf = [it / self.VEL for it in vf]
        
        # Assembling the constraints
        ceq = [(self.ta.state[0] - rf[0])]
        ceq += [(self.ta.state[1] - rf[1])]
        ceq += [(self.ta.state[2] - rf[2])]
        ceq += [(self.ta.state[3] - vf[0])]
        ceq += [(self.ta.state[4] - vf[1])]
        ceq += [(self.ta.state[5] - vf[2])]
        ceq += [self.ta.state[13]]  # lm = 0
        ceq += [sum([it * it for it in x[:8]]) - 1.0]  # |lambdas|^2 = 1
        return [1.0] + ceq

    def gradient(self, x):
        # Single Shooting of variational equations
        self.set_ta_var_state(x[:8])
        self.ta_var.propagate_until(x[8])
        # Computing the target position, velocity and acceleration at epoch
        # NOTE: if the planet is not keplerian the acceleration will be assumed as Keplerian
        rf, vf = self.target.eph(self.t0 + x[8] * self.TIME * _pk.SEC2DAY)
        vf = [it / self.VEL for it in vf]
        rf = [it / self.L for it in rf]
        af = -self.mu / _np.linalg.norm(rf) ** 3 * _np.array(rf)
        # Assembling the gradient
        # Fitness (sparsity takes care of the fitness since it does not depend on the dv)
        retval = []
        dyn = self.dyn_func(self.ta_var.state[:14], pars=self.ta_var.pars[:2])
        # Constraints
        for i in range(3):
            sl = self.ta_var.get_vslice(order=1, component=i)
            retval.extend(list(self.ta_var.state[sl]))
            retval.append(dyn[i] - vf[i])
        for i in range(3, 6):
            sl = self.ta_var.get_vslice(order=1, component=i)
            retval.extend(list(self.ta_var.state[sl]))
            retval.append(dyn[i] - af[i - 3])
        # Constraint on mass costate
        sl = self.ta_var.get_vslice(order=1, component=13)
        retval.extend(list(self.ta_var.state[sl]))
        retval.append(dyn[6])
        # Norm constraint
        for item in x[:-1]:
            retval.extend([2 * item])
        return retval

    def gradient_sparsity(self):
        # x = [lx,ly,lz,lvx,lvy,lvz,lm,l0,tof]
        retval = []
        for i in range(1, 8):
            for j in range(9):
                retval.append((i, j))
        # The norm constraint does not depend on tof
        for i in range(8, 9):
            for j in range(8):
                retval.append((i, j))
        return retval

    def has_gradient(self):
        return self.with_gradient

    def get_nec(self):
        return 8

    def plot_trajectory(self, x, N=100, ax3D=None):
        """
        This function plots the trajectory encoded in the decision vector x.

        Args:
            *x* (:class:`list`): the decision vector.
            *N* (:class:`int`): the number of points to use in the plot.
            *ax3D* (:class:`matplotlib.axes._axes.Axes`): the axis to use for the plot. If None, a new axis is created.

        Returns:
            *ax3D* (:class:`matplotlib.axes._axes.Axes`): the axis of the plot.
        """
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d.art3d import Line3DCollection
        import matplotlib.colors as mcolors

        # Single Shooting
        self.set_ta_state(x[:8])
        t_grid = _np.linspace(0, x[8], N)
        sol = self.ta.propagate_grid(t_grid)
        # We make the axis if needed
        if ax3D is None:
            ax3D = _pk.plot.make_3Daxis()
        # Adding the main body
        _pk.plot.add_sun(ax3D)
        _pk.plot.add_planet_orbit(ax3D, self.source, c="gray", units=_pk.AU)
        _pk.plot.add_planet_orbit(ax3D, self.target, c="gray", units=_pk.AU)
        # Plotting the trajectory
        # Assuming sol[-1] is an array of shape (N, 3)
        state = sol[-1]  # Extract (x, y, z) coordinates
        xx, yy, zz = state[:, 0], state[:, 1], state[:, 2]

        ax3D.plot(xx, yy, zz)

        # Plotting the boundary conditions
        ax3D.scatter(self.posvel0[0][0], self.posvel0[0][1], self.posvel0[0][2])
        rf, _ = self.target.eph(self.t0 + x[8] * self.TIME * _pk.SEC2DAY)
        rf = [it / self.L for it in rf]
        ax3D.scatter(rf[0], rf[1], rf[2])
        return ax3D

    def plot_misc(self, x, N=100, **kwargs):
        """
        This function plots the throttle, switching function and Hamiltonian
        of the trajectory encoded in the decision vector x.

        Args:
            *x* (:class:`list`): the decision vector.
            *N* (:class:`int`): the number of points to use in the plot.
            *\\*\\*kwargs*: Additional keyword arguments to pass to the plt.subplots function.

        Returns:
            *axs* (:class:`list`): the list of axis of the plots.
        """
        from matplotlib import pyplot as plt

        # Single Shooting
        self.set_ta_state(x[:8])
        t_grid = _np.linspace(0, x[8], N)
        sol = self.ta.propagate_grid(t_grid)
        # Retreive useful cfuncs
        # The Hamiltonian
        H_func = _pk.ta.get_pc_H_cfunc(_pk.optimality_type.TIME)
        # The switching function
        SF_func = _pk.ta.get_pc_SF_cfunc(_pk.optimality_type.TIME)
        # The thrust direction
        i_vers_func = _pk.ta.get_pc_i_vers_cfunc(_pk.optimality_type.TIME)

        # Create axis
        _, axs = plt.subplots(2, 2, **kwargs)
        axs[0, 0].set_title("Mass")
        axs[0, 0].plot(t_grid, self.MASS * sol[-1][:, 6].T)
        # Plot thrust direction
        thrust_dir = i_vers_func(_np.ascontiguousarray(sol[-1][:, 10:13].T))

        for i in range(3):
            axs[0, 1].plot(t_grid, thrust_dir[i, :])
        axs[0, 1].set_title("Thrust direction")

        # Plot Hamiltonian
        all_pars = list(self.ta.pars) + [1., x[-2]]
        Ham = H_func(
            _np.ascontiguousarray(sol[-1].T),
           pars=_np.ascontiguousarray(_np.tile(all_pars, (N, 1)).T),
         )
        axs[1, 0].set_title("Hamiltonian")
        axs[1, 0].plot(t_grid, _np.squeeze(Ham))
        
        # Plot switching function
        SF = SF_func(
            _np.ascontiguousarray(sol[-1].T),
            pars=_np.ascontiguousarray(_np.tile(all_pars, (N, 1)).T),
        )
        axs[1, 1].set_title("Switching Function")
        axs[1, 1].plot(t_grid, _np.squeeze(SF))
        axs[1, 1].hlines(0, t_grid[0], t_grid[-1], "k")

        plt.tight_layout()

        return axs


del _posvel0, _posvelf, _pl0, _plf

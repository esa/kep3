import numpy as _np
import heyoka as _hy

class zoh:
    def __init__(
        self,
        state0,
        controls,
        state1,
        tgrid,
        cut,
        tas,
    ):
        # We store the constructor args
        self.state0 = state0
        self.controls = controls
        self.state1 = state1
        self.tgrid = tgrid
        self.cut = cut
        # Store the tas
        self.ta = tas[0]
        self.ta_var = tas[1]
        # And save the non control parameter values for cfunc calls
        self.pars_no_control = self.ta.pars[4:].tolist()

        # We compute convenient quantities
        self.nseg = len(self.controls) // 4
        self.nseg_fwd = int(self.nseg * cut)
        self.nseg_bck = self.nseg - self.nseg_fwd

        # Sanity checks
        # On the tas
        if len(self.ta.state) != 7:
            raise ValueError(
                f"Attempting to construct a zoh_leg with a Taylor Adaptive integrator state dimension of {len(self.ta.state)}, while 7 is required"
            )
        if len(self.ta.pars) <= 4:
            raise ValueError(
                f"Attempting to construct a zoh_leg with a Taylor Adaptive integrator parameters dimension of {len(self.ta.pars)}, while >=4 is required"
            )
        if len(self.ta_var.state) != 7 + 7 * 7 + 7 * 4:
            raise ValueError(
                f"Attempting to construct a zoh_leg with a variational Taylor Adaptive integrator state dimension of {len(self.ta_var.state)}, while 84 is required"
            )
        if len(self.ta_var.pars) != len(self.ta.pars):
            raise ValueError(
                f"While constructing a zoh_leg, the number of parameters in the variational version of the Taylor integrator and the non variational version were detected as different, while its required they are equal"
            )

        # On the rest
        if len(self.controls) % 4 > 0:
            raise ValueError(
                "In a zoh_leg controls must be multiple of 4: [T, ix, iy, iz] * nseg"
            )
        if len(tgrid) != self.nseg + 1:
            raise ValueError(
                "The t_grid and the controls have incompatible lenghts. It must be nseg*4 and nseg+1"
            )

        # this assumes state of 7 and 4 pars
        self.ic_var = (_np.hstack((_np.eye(7, 7), _np.zeros((7, 4))))).flatten()

        # We compile the function for the dynamics (this is used in the gradient computations)
        sys = self.ta.sys
        vars = [it[0] for it in sys]
        dyn = [it[1] for it in sys]
        self.dyn_cfunc = _hy.cfunc_dbl(dyn, vars, compact_mode=True)

    def compute_mismatch_constraints(self):
        # Forward segments
        self.ta.time = self.tgrid[0]
        self.ta.state[:] = self.state0
        for i in range(self.nseg_fwd):
            # setting T, ix, iy, iz
            self.ta.pars[0:4] = self.controls[4 * i : 4 * i + 4]
            # propagating
            self.ta.propagate_until(self.tgrid[i + 1])
        state_fwd = self.ta.state.copy()

        # Backward segments
        self.ta.time = self.tgrid[-1]
        self.ta.state[:] = self.state1
        for i in range(self.nseg_bck):
            # setting T, ix, iy, iz
            self.ta.pars[0] = self.controls[-4 * i - 4]
            self.ta.pars[1] = self.controls[-4 * i - 3]
            self.ta.pars[2] = self.controls[-4 * i - 2]
            self.ta.pars[3] = self.controls[-4 * i - 1]
            # propagating
            self.ta.propagate_until(self.tgrid[-2 - i])
        state_bck = self.ta.state
        return (state_fwd - state_bck).tolist()

    def compute_throttle_constraints(self):
        retval = [0] * self.nseg
        for i in range(self.nseg):
            retval[i] = (
                self.controls[1 + 4 * i] ** 2
                + self.controls[2 + 4 * i] ** 2
                + self.controls[3 + 4 * i] ** 2
                - 1
            )
        return retval

    def compute_tc_grad(self):
        retval = _np.zeros((self.nseg, self.nseg * 4))
        for i in range(self.nseg):
            retval[i, 4 * i + 1] = 2 * self.controls[4 * i + 1]
            retval[i, 4 * i + 2] = 2 * self.controls[4 * i + 2]
            retval[i, 4 * i + 3] = 2 * self.controls[4 * i + 3]
        return retval

    def compute_mc_grad(self):
        # We will return 4 matrices: dmc/dx0, dmc/dx1, dmc/dcontrol, dmc/dtgrid
        # STMs -> forward
        M_seg_fwd = []  # M10, M21, M32, ....
        M_fwd = []  # Mf0, Mf1, Mf2, ...
        C_seg_fwd = []  # dx_(i+1)/dcontrols_i
        C_fwd = _np.zeros((7, 4 * self.nseg_fwd))  # dx_f/dcontrols_i
        dyn_fwd = []
        self.ta_var.time = self.tgrid[0]
        self.ta_var.state[:7] = self.state0
        for i in range(self.nseg_fwd):
            self.ta_var.state[7:] = self.ic_var
            # setting T, ix, iy, iz
            self.ta_var.pars[0:4] = self.controls[4 * i : 4 * i + 4]
            # propagating
            self.ta_var.propagate_until(self.tgrid[i + 1])
            # extracting the segment STMs
            M_seg_fwd.append(self.ta_var.state[7:].reshape(7, 11)[:, :7].copy())  # 7x7
            # extracting the control sensitivities in across the single segment
            C_seg_fwd.append(self.ta_var.state[7:].reshape(7, 11)[:, 7:].copy())  # 7x4
            # computing the dynamics for future usage
            dyn_fwd.append(
                self.dyn_cfunc(
                    self.ta_var.state[:7],
                    pars=self.controls[4 * i : 4 * i + 4] + self.pars_no_control,
                )
            )
        # We compute the STMs - Mf0, Mf1, Mf2, ...
        cur = _np.eye(7)
        for M in reversed(M_seg_fwd):
            cur = cur @ M
            M_fwd.append(cur)
        M_fwd = list(reversed(M_fwd)) + [_np.eye(7)]

        # 1 - dmc/dx0
        dmc_dx0 = M_fwd[0]

        # 2 - dmc/dcontrols
        i = 0
        for M, C in zip(M_fwd[1:], C_seg_fwd):
            C_fwd[:, 4 * i : 4 * i + 4] = M @ C
            i += 1

        # 3 - dmc/dtgrid
        dmcdtgrid = _np.zeros((7, self.nseg + 1))
        if self.nseg_fwd > 0:
            dmcdtgrid[:, 0] = -M_fwd[1] @ dyn_fwd[0]
            dmcdtgrid[:, i] = M_fwd[-1] @ dyn_fwd[-1]

            for i in range(1, self.nseg_fwd):
                dmcdtgrid[:, i] = M_fwd[i + 1] @ (
                    M_seg_fwd[i] @ dyn_fwd[i - 1] - dyn_fwd[i]
                )

        # STMs -> backward
        M_seg_bck = []  # M10, M21, M32, .... note: numbering reversed
        M_bck = []  # Mf0, Mf1, Mf2, .... note: numbering reversed
        C_seg_bck = []  # dxi/dcontrols
        C_bck = _np.zeros((7, 4 * self.nseg_bck))
        dyn_bck = []
        self.ta_var.time = self.tgrid[-1]
        self.ta_var.state[:7] = self.state1
        for i in range(self.nseg_bck):
            self.ta_var.state[7:] = self.ic_var
            # setting T, ix, iy, iz
            self.ta_var.pars[0] = self.controls[-4 * i - 4]
            self.ta_var.pars[1] = self.controls[-4 * i - 3]
            self.ta_var.pars[2] = self.controls[-4 * i - 2]
            self.ta_var.pars[3] = self.controls[-4 * i - 1]
            # propagating
            self.ta_var.propagate_until(self.tgrid[-2 - i])
            # extracting the segment STMs
            M_seg_bck.append(self.ta_var.state[7:].reshape(7, 11)[:, :7].copy())
            # extracting the control sensitivities in across the single segment
            C_seg_bck.append(self.ta_var.state[7:].reshape(7, 11)[:, 7:].copy())  # 7x4
            # computing the dynamics for future usage
            dyn_bck.append(
                self.dyn_cfunc(
                    self.ta_var.state[:7],
                    pars=self.controls[-4 * i - 4 : -4 * i - 1]
                    + [self.controls[-4 * i - 1]]
                    + self.pars_no_control,
                )
            )

        # We compute the STMs - Mf0, Mf1, Mf2, ...
        cur = _np.eye(7)
        for M in reversed(M_seg_bck):
            cur = cur @ M
            M_bck.append(cur)
        M_bck = list(reversed(M_bck)) + [_np.eye(7)]

        # 1 - dmc/dxf
        dmc_dx1 = -M_bck[0]  # mc = xfwd-x_bck hence the minus

        # 2 - dmc/dcontrols
        # NOTE: the controls are reversed in order and we need to account for it
        i = 0
        for M, C in zip(M_bck[1:], C_seg_bck):
            C_bck[:, 4 * self.nseg_bck - 4 - 4 * i : 4 * self.nseg_bck - 4 * i] = M @ C
            i += 1

        # 3 - dmc/dtgrid
        if (
            self.nseg_bck > 0
        ):  # we skip this in the corner case where no bck seg are there
            dmcdtgrid[:, -1] = M_bck[1] @ dyn_bck[0]
            dmcdtgrid[:, self.nseg_fwd] -= (
                M_bck[-1] @ dyn_bck[-1]
            )  # This is for the mid time point shared fwd and bck: we need to subtract to the existing

            for i in range(1, self.nseg_bck):
                dmcdtgrid[:, -1 - i] = -M_bck[i + 1] @ (
                    M_seg_bck[i] @ dyn_bck[i - 1] - dyn_bck[i]
                )

        # dmc/dcontrols (7x4*nseg)
        dmc_dcontrols = _np.hstack((C_fwd, -C_bck))  # mc = xfwd-x_bck hence the minus

        return dmc_dx0, dmc_dx1, dmc_dcontrols, dmcdtgrid

    def add_to_3Daxis(self, ax, N=100, scale=1.0, plot_time_grid=False):
        # Forward segments
        self.ta.time = self.tgrid[0]
        self.ta.state[:] = self.state0
        for i in range(self.nseg_fwd):
            # setting T, ix, iy, iz
            self.ta.pars[0:4] = self.controls[4 * i : 4 * i + 4]
            # propagating
            plot_grid_fwd = _np.linspace(self.tgrid[i], self.tgrid[i + 1], N)
            sol_fwd = self.ta.propagate_grid(plot_grid_fwd)[-1]
            ax.plot(
                sol_fwd[:, 0] * scale,
                sol_fwd[:, 1] * scale,
                sol_fwd[:, 2] * scale,
                c="blue",
            )
            if plot_time_grid:
                ax.scatter(
                    sol_fwd[-1, 0] * scale,
                    sol_fwd[-1, 1] * scale,
                    sol_fwd[-1, 2] * scale,
                    c="blue",
                )

        # Backward segments
        self.ta.time = self.tgrid[-1]
        self.ta.state[:] = self.state1
        for i in range(self.nseg_bck):
            # setting T, ix, iy, iz
            self.ta.pars[0] = self.controls[-4 * i - 4]
            self.ta.pars[1] = self.controls[-4 * i - 3]
            self.ta.pars[2] = self.controls[-4 * i - 2]
            self.ta.pars[3] = self.controls[-4 * i - 1]
            # propagating
            plot_grid_bck = _np.linspace(self.tgrid[-1 - i], self.tgrid[-2 - i], N)
            sol_bck = self.ta.propagate_grid(plot_grid_bck)[-1]
            ax.plot(
                sol_bck[:, 0] * scale,
                sol_bck[:, 1] * scale,
                sol_bck[:, 2] * scale,
                c="darkorange",
            )
            if plot_time_grid:
                ax.scatter(
                    sol_bck[-1, 0] * scale,
                    sol_bck[-1, 1] * scale,
                    sol_bck[-1, 2] * scale,
                    c="darkorange",
                )

        # Adding start and final position
        ax.scatter(
            self.state0[0] * scale,
            self.state0[1] * scale,
            self.state0[2] * scale,
            c="blue",
        )
        ax.scatter(
            self.state1[0] * scale,
            self.state1[1] * scale,
            self.state1[2] * scale,
            c="darkorange",
        )
        return ax
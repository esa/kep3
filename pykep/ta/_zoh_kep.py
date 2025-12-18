import heyoka as _hy
import numpy as _np
from copy import deepcopy as _deepcopy


def zoh_kep_dyn():
    """
    The dynamics in Cartesian coordinates of a constant thrust mass-varying spacecraft
    orbiting a main central body.
    
    We consider the motion of a spacecraft of mass :math:`m` with a position 
    :math:`\\mathbf{r}` and velocity :math:`\\mathbf{v}` subject only to a main body
    gravitational attraction in some inertial reference frame. The spacecraft also
    has an ion thruster characterized by :math:`v_{eff} = I_{sp}g_0` and thrusting as
    :math:`\mathbf T = T [i_x, i_y, i_z]`. 
    
    The dynamics are thus given by:
    
    .. math::
        \\begin{array}{l}
        \\dot{\\mathbf r} = \\mathbf v \\\\
        \\dot{\\mathbf v} = -\\frac{1}{r^3} \\mathbf r + \\frac{T}{m} \\mathbf {\\hat i} \\\\
        \\dot m = - c T 
        \\end{array}
        
    where, :math:`c = \frac 1{v_{eff}}`. 
    The gravitational parameter is assumed as unitary, all other quantities to be 
    scaled accordingly.
    
    State in heyoka is :math:`[x,y,z,vx,vy,vz,m]`
    Parameters in heyoka are: :math:`[T, ix, iy, iz] + [c]`
    
    Returns:
        :class:`list` [ :class:`tuple` (:class:`hy::expression`, :class:`hy::expression` )]: The dynamics in the form [(p, dp), ...]

    """
    # Estalishing the state
    x, y, z, vx, vy, vz, m = _hy.make_vars("x", "y", "z", "vx", "vy", "vz", "m")

    # Naming the system parametes
    c = _hy.par[4]

    # Naming the system controls
    T = _hy.par[0]
    ix, iy, iz = _hy.par[1], _hy.par[2], _hy.par[3]

    # Auxiliary expressions
    r2 = _hy.sum([x**2, y**2, z**2])

    # Building the dynamics
    xdot = vx
    ydot = vy
    zdot = vz
    vxdot = -1.0 * pow(r2, -1.5) * x + T * ix / m
    vydot = -1.0 * pow(r2, -1.5) * y + T * iy / m
    vzdot = -1.0 * pow(r2, -1.5) * z + T * iz / m
    mdot = -c * T
    retval = [
        (x, xdot),
        (y, ydot),
        (z, zdot),
        (vx, vxdot),
        (vy, vydot),
        (vz, vzdot),
        (m, mdot),
    ]
    return retval


# We mimick in python the C++ global caching mechanism for taylor_adaptive
# instances, so that we do not create multiple instances with
# the same tolerance.
_ta_zoh_kep_cache = dict()


def get_zoh_kep(tol: float):
    """
    Returns a Taylor adaptive propagator (Heyoka) for the zoh_kep_dyn dynamics retreiving
    one from a global cache.

    This is the initial value problem of a constant thrust mass-varying spacecraft orbiting a primary.
    Thrust direction is fixed in the inertial frame.

    If the requested propagator was never created this will create it, else it will
    return the one from the global cache (not a copy, user will need to copy it
    if needed).

    The dynamics is that returned by :func:`~pykep.ta.zoh_kep_dyn`.

    Args:
        *tol* (:class:`float`): the tolerance of the Taylor adaptive propagator.

    Returns:
        :class:`hy::taylor_adaptive`: The Taylor adaptive propagator.

    Examples:
    >>> import pykep as pk
    >>> ta = pk.ta.zoh_kep_dyn(tol = 1e-16)
    >>> ta.time = 0.
    >>> ta.state[:] = [1.2,0.1,0.1,0.1,1.,0.123,1.]
    >>> veff = 1.32
    >>> controls = [0.022, 0.023, -0.21, 0.1]
    >>> tof = 1.23
    >>> ta.pars[:] = controls + [1 / veff]
    >>> ta.propagate_until(tof)
    """
    # Lookup.
    if _ta_zoh_kep_cache.get(tol) is None:
        # Cache miss, create new one.
        init_state = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        new_ta = _hy.taylor_adaptive(
            zoh_kep_dyn(), init_state, tol=tol, pars=[1.0, 1.0, 0.0, 0.0, 0.0]
        )
        _ta_zoh_kep_cache[tol] = new_ta
        return new_ta
    else:
        # Cache hit, return existing.
        return _ta_zoh_kep_cache[tol]


_ta_zoh_kep_var_cache = dict()


def get_zoh_kep_var(tol: float):
    """Returns a (order 1) variational Taylor adaptive propagator (Heyoka)
    for the zoh_kep dynamics retreiving one from a global cache if possible.

    .. note:
    Variations are only considered with repsect to initial conditions and the
    thrust parameters :math:`[T, i_x, i_y, i_z]`.

    The dynamics is that returned as expressions by :func:`~pykep.ta.zoh_kep_dyn`
    and also used in the non variational version returend by :func:`~pykep.ta.get_zoh_kep`

    Args:
        *tol* (:class:`float`): the tolerance of the variational Taylor adaptive propagator.

    Returns:
        :class:`hy::taylor_adaptive`: The variational Taylor adaptive propagator.

    Examples:
    >>> import pykep as pk
    >>> veff = 1.32
    >>> controls = [0.022, 0.023, -0.21, 0.1]
    >>> tof = 1.23
    >>> ta = pk.ta.get_zoh_kep_var(tol = 1e-16)
    >>> ta.time = 0.
    >>> ta.state[:7] = [1.2,0.1,0.1,0.1,1.,0.123,1.]
    >>> ta.pars[:] = controls + [1 / veff]
    >>> ta.propagate_until(tof)
    """
    # Lookup.
    if _ta_zoh_kep_var_cache.get(tol) is None:
        # Cache miss, create new one.
        x, y, z, vx, vy, vz, m = _hy.make_vars("x", "y", "z", "vx", "vy", "vz", "m")
        vsys = _hy.var_ode_sys(
            zoh_kep_dyn(),
            [x, y, z, vx, vy, vz, m, _hy.par[0], _hy.par[1], _hy.par[2], _hy.par[3]],
            1,
        )
        new_ta = _hy.taylor_adaptive(
            vsys,
            tol=tol,
            compact_mode=True,
        )
        # Insert in cache and return
        _ta_zoh_kep_var_cache[tol] = new_ta
        return new_ta
    else:
        # Cache hit, return existing.
        return _ta_zoh_kep_var_cache[tol]

from copy import deepcopy as _deepcopy
from math import log, exp, cos, acos, sin, pi, sqrt, atan2, asin

def alpha2direct(alphas, T):
    """alpha2direct(x)

    Args:
        *alphas* (``array-like``): a sequence of transfer times encoded using the alpha encoding.

        *T* (:class:`float`): the total transfer time.

    Returns:
        :class:`list`:: The encoded transfer times
    """
    log_alphas = [log(item) for item in alphas]
    retval = [it / sum(log_alphas) * T for it in log_alphas]
    return retval


def direct2alpha(x):
    """direct2alpha(x)

    Args:
        *x* (``array-like``): a sequence of transfer times.

    Returns:
        :class:`list`:, :class:`float`: The alpha-encoded transfer times, the total transfer time (for cenvenience)
    """
    T = sum(x)
    retval = [exp(it / (-T)) for it in x]
    return retval, T


def eta2direct(x, max_tof):
    """eta2direct(x)

    Args:
        *x* (``array-like``): a sequence of transfer times encoded using the eta encoding.

        *max_tof* (:class:`float`): maximum time of flight (might actually be less according to the value of the last eta.)

    Returns:
        :class:`list`: The encoded transfer times
    """
    n = len(x)
    # we assemble the times of flight
    T = [0] * n
    T[0] = max_tof * x[0]
    for i in range(1, len(T)):
        T[i] = (max_tof - sum(T[:i])) * x[i]
    return T


def direct2eta(x, max_tof):
    """direct2eta(x)

    Args:
        *x* (``array-like``):  a sequence of transfer times

        *max_tof* (:class:`float`): maximum time of flight (might actually be less according to the value of the last eta.)

    Returns:
        :class:`list`: The eta-encoded transfer times
    """
    retval = _deepcopy(x)
    retval[0] = x[0] / max_tof
    for i in range(1, len(x)):
        retval[i] = x[i] / (max_tof - sum(x[:i]))
    return retval


def uvV2cartesian(uvV):
    """This function converts the uvV encoding of a vector to cartesian coordinates

    Args:
        *uvV* (``array-like``): a sequence of 3 floats representing the vector in uvV encoding.

    Returns:
        :class:`list`:: The vector in cartesian coordinates
    """
    u, v, V = uvV
    theta = 2 * pi * u
    tmp = 2 * v - 1
    # Protecting against nans
    if (tmp*tmp) > 1:
        if tmp < -1:
            phi = pi/2
        elif tmp > 1:
            phi = -pi/2
    else:
        phi = acos(2 * v - 1) - pi / 2

    Vinfx = V * cos(phi) * cos(theta)
    Vinfy = V * cos(phi) * sin(theta)
    Vinfz = V * sin(phi)
    return [Vinfx, Vinfy, Vinfz]


def cartesian2uvV(V):
    """This function converts the cartesian coordinates of a vector to uvV encoding

    Args:
        *V* (``array-like``): a sequence of 3 floats representing the vector in cartesian coordinates.

    Returns:
        :class:`list`:: The vector in uvV encoding
    """
    v_norm = sqrt(V[0] ** 2 + V[1] ** 2 + V[2] ** 2)
    theta = atan2(V[1], V[0])
    sin_phi = V[2] / v_norm
    return [theta / 2 / pi, (-sin_phi + 1) / 2, v_norm]

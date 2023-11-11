import pykep as _pk
import numpy as _np

def add_planet(ax, pla: _pk.planet, when, label=None, c="gray", s=10, units=_pk.AU):
    """Adds a planet to *ax*.

    Args:
        *ax* (:class:`mpl_toolkits.mplot3d.axes3d.Axes3D`): the 3D axis.

        *pla* (:class:`~pykep.planet`): the planet.

        *when* (:class:`~pykep.epoch` or :class:`float`): the epoch (in mjd2000 if float).

        *label* (:class:`str`, optional): the plot label. Defaults to the planet name as retuned by the corresponding planet method.

        *c* (:class:`str`, optional): the color. Defaults to "gray".

        *s* (:class:`int`, optional): the size. Defaults to 10.

        *units* (:class:`float`, optional): length units to be used. Defaults to pk.AU.

    Returns:
        :class:`mpl_toolkits.mplot3d.axes3d.Axes3D`: the 3D axis.

    Examples:
        >>> import pykep as pk
        >>> ax = pk.plot.make_3Daxis()
        >>> pk.plot.add_sun(ax, label="Sun")
        >>> udpla = pk.udpla.jpl_lp(body="EARTH")
        >>> earth = pk.planet(udpla)
        >>> pk.plot.add_planet_orbit(ax, earth, plot_range = [0, 365.25], c = "royalblue", label = "")
        >>> pk.plot.add_planet(ax, earth, when = pk.epoch(0), c = "royalblue")
    """
    # When no label is passed, the default behaviour is to add one with the planet name as returned by the corresponding method.
    if label is None:
        label = pla.get_name()

    # We compute the ephemerides
    r, _ = pla.eph(when)

    # And plot them as a scatter3D
    ax.scatter(r[0] / units, r[1] / units, r[2] / units, s=s, c=c, label=label)

    # Returning the axes.
    return ax


def add_planet_orbit(
    ax,
    pla: _pk.planet,
    plot_range=None,
    label=None,
    c="gray",
    units=_pk.AU,
    N=60,
):
    """Adds a planet orbit to *ax*.

    Args:
        *ax* (:class:`mpl_toolkits.mplot3d.axes3d.Axes3D`): the 3D axis.

        *pla* (:class:`~pykep.planet`): the planet.

        *plot_range* (:class:`list`, optional): the starting and end mjd2000 to be plotted. Defaults to [0., one_orbital_period] if a period can be computed.

        *label* (:class:`str`, optional): the plot label. Defaults to the planet name as retuned by the corresponding planet method.

        *c* (:class:`str`, optional): the color. Defaults to "gray".

        *s* (:class:`int`, optional): the size. Defaults to 10.

        *units* (:class:`float`, optional): length units to be used. Defaults to pk.AU.

    Returns:
        :class:`mpl_toolkits.mplot3d.axes3d.Axes3D`: the 3D axis.

    Examples:
        >>> import pykep as pk
        >>> ax = pk.plot.make_3Daxis()
        >>> pk.plot.add_sun(ax, label="Sun")
        >>> udpla = pk.udpla.jpl_lp(body="EARTH")
        >>> earth = pk.planet(udpla)
        >>> pk.plot.add_planet_orbit(ax, earth, plot_range = [0, 365.25], c = "royalblue", label = "")
        >>> pk.plot.add_planet(ax, earth, when = pk.epoch(0), c = "royalblue")
    """
    # When no label is passed, the default behaviour is to add one with the planet name as returned by the corresponding method.
    if label is None:
        label = pla.get_name()

    # If the plot_range is not defined by the user, then a defult is attempted [0,T]
    if plot_range is None:
        try:
            T = pla.period() * _pk.SEC2DAY
        except NotImplementedError as e:
            print(
                f"PyKEP ERROR: Cannot compute the orbital period when plotting the planet {pla.get_name()}, please define an explicit plot_range."
            )
            raise e
        epochs = _np.linspace(0.0, T, N)
    else:
        epochs = _np.linspace(plot_range[0], plot_range[1], N)

    # Compute the ephemerides
    rvs = pla.eph_v(epochs)[:, :3] / units

    # Plots the planet in the range.
    ax.plot(rvs[:, 0], rvs[:, 1], rvs[:, 2], c=c, label=label)

    # Returning the axes.
    return ax

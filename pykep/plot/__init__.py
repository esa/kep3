# We import symbols explicitly with an underscore to hide them from the
# imported symbols
import matplotlib.pyplot as _plt
from mpl_toolkits.mplot3d import axes3d as _axes3d

from ._planet import add_planet_orbit, add_planet, add_solar_system

def make_3Daxis(**kwargs):
    """Constructs and returns a 3D axis.  All kwargs are forwarded to the call to `figure()` in matplotlib.

    Returns:
        :class:`mpl_toolkits.mplot3d.axes3d.Axes3D`: the 3D axis
    """
    ax = _plt.figure(**kwargs).add_subplot(projection='3d')
    return ax

def add_sun(ax, **kwargs):
    """Adds the Sun to *ax*. All kwargs are forwarded to the scatter method of *ax*.

    Args:
        *ax* (:class:`mpl_toolkits.mplot3d.axes3d.Axes3D`): the 3D axis.

    Returns:
        :class:`mpl_toolkits.mplot3d.axes3d.Axes3D`: the 3D axis.
    """
    kwargs.setdefault('c', 'y')
    kwargs.setdefault('s', 30)
    kwargs.setdefault('label', 'Sun')

    ax.scatter(0,0,0, **kwargs)
    return ax



import pykep as pk
import numpy as np

# Earth
udpla = pk.udpla.jpl_lp(body="EARTH")
earth = pk.planet(udpla)

sf = pk.leg.sims_flanagan()
nseg = 20
sf.throttles=[0,0,0] * nseg
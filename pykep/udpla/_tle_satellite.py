from sgp4.api import Satrec
from sgp4 import exporter
import pykep as pk

class tle_satellite:
    """__init__(line1, line2)

    This User Defined Planet (UDPLA) represents a satellite orbiting the Earth and defined in the TLE format
    and propagated using the SGP4 propagator.

    .. note::
       The resulting ephemerides will be returned in SI units and in the True Equator Mean Equinox (TEME) reference frame

    """
    def __init__(self, line1, line2):
        import pykep as pk
        self.satellite = Satrec.twoline2rv(line1, line2)
        self.e = 0
        self.ref_epoch = pk.epoch(self.satellite.jdsatepoch + self.satellite.jdsatepochF, pk.epoch.julian_type.JD)
    def eph(self, ep):
        jd = ep.jd()
        jd_i = int(jd)
        jd_fr = jd-jd_i
        self.e, r, v = self.satellite.sgp4(jd_i, jd_fr)
        return [[it*1000 for it in r], [it*1000 for it in v]]
    def get_name(self):
        return self.satellite.satnum_str + " - SGP4"
    def get_extra_info(self):
        line1, line2 = exporter.export_tle(self.satellite)
        return "TLE line1: " + line1 + "\nTLE line2: " + line2 
    def get_mu_central_body(self):
        return pk.MU_EARTH
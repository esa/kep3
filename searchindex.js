Search.setIndex({"docnames": ["anomalies", "api", "constants", "elements", "epoch", "index", "lambert", "notebooks/anomalies", "notebooks/benchmarkudplaiface", "notebooks/epochs", "notebooks/interface_to_spice", "notebooks/planet", "notebooks/propagate_lagrangian", "planet", "propagation", "tutorials", "udplas"], "filenames": ["anomalies.rst", "api.rst", "constants.rst", "elements.rst", "epoch.rst", "index.md", "lambert.rst", "notebooks/anomalies.ipynb", "notebooks/benchmarkudplaiface.ipynb", "notebooks/epochs.ipynb", "notebooks/interface_to_spice.ipynb", "notebooks/planet.ipynb", "notebooks/propagate_lagrangian.ipynb", "planet.rst", "propagation.rst", "tutorials.rst", "udplas.rst"], "titles": ["Anomalies Conversions", "API", "Global constants", "Orbital Elements", "Epoch class", "Welcome to pykep\u2019s documentation!", "Lambert class", "The various anomalies in pykep", "&lt;no title&gt;", "Epochs and Julian Dates", "Interfacing to SPICE and JPL DE ephs", "Ephemerides", "Lagrange Propagation", "Planet class", "Numerical Propagation", "Tutorials", "List of user implemented planets (UDPLAs)"], "terms": {"In": [0, 2, 3, 4, 9, 10, 11, 12, 13], "pykep": [0, 1, 3, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16], "we": [0, 4, 7, 9, 10, 11, 12, 13], "adopt": 0, "follow": [0, 9, 10, 11, 12, 13], "name": [0, 2, 7, 10, 11, 13, 16], "variou": [0, 2, 3, 11, 15], "m": [0, 2, 3, 5, 7, 10, 12, 13], "i": [0, 1, 3, 4, 5, 7, 9, 10, 11, 12, 13, 14, 16], "mean": [0, 3, 7, 11, 16], "e": [0, 3, 7, 8, 9, 13, 16], "eccentr": [0, 7, 11, 12], "l": [0, 3, 7, 12], "true": [0, 3, 6, 7, 9, 11, 12, 13, 14, 16], "longitud": [0, 3, 7], "lambda": [0, 7], "h": [0, 3, 7], "hyperbol": [0, 7, 12], "n": [0, 2, 7, 12], "zeta": [0, 7], "gudermannian": [0, 7], "function": [0, 2, 3, 7, 9, 10, 12, 13, 16], "variabl": [0, 6, 7], "symbol": [0, 2, 7], "can": [0, 2, 4, 7, 9, 10, 11, 12, 13, 15, 16], "spell": [0, 12], "out": [0, 12], "like": [0, 6, 14, 16], "made": [0, 6, 10, 11], "lowercas": 0, "below": [0, 12], "list": [0, 1, 5, 10, 13, 14], "allow": [0, 5, 6, 9, 11, 12, 13], "convert": [0, 3, 7, 9, 10], "from": [0, 3, 4, 7, 8, 9, 10, 12, 13, 16], "one": [0, 2, 4, 9, 10, 12, 13], "anoth": 0, "version": [0, 3, 7, 12, 13], "m2e": [0, 7], "ecc": [0, 7], "requir": [0, 10], "1": [0, 2, 3, 4, 6, 7, 8, 10, 11, 12, 13, 14, 16], "arg": [0, 4, 6, 13, 14, 16], "float": [0, 4, 6, 9, 13, 14, 16], "rad": 0, "return": [0, 4, 11, 13, 14, 16], "pi": [0, 3, 6, 7, 12, 13, 14], "exampl": [0, 4, 6, 7, 12, 13, 14, 16], "import": [0, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16], "pk": [0, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16], "2": [0, 2, 3, 6, 7, 8, 10, 12, 13, 14, 16], "0": [0, 2, 3, 5, 6, 7, 8, 9, 11, 12, 13, 14, 16], "296254963787226": 0, "e2m": 0, "5": [0, 7, 8, 9, 12], "4520574461395797": 0, "m2f": 0, "32": [0, 9], "65": 0, "4497431281728277": 0, "f2m": 0, "f": [0, 3, 7, 8, 10, 12, 16], "34": 0, "67": 0, "05065883735669101": 0, "e2f": 0, "5502639747136633": 0, "f2e": 0, "1082931139529482": 0, "n2h": [0, 7], "10": [0, 2, 9, 10, 16], "12836469743916526": 0, "h2n": 0, "14": [0, 16], "377641187853621": 0, "n2f": 0, "13": [0, 4], "45": [0, 12], "7373697968359353": 0, "f2n": 0, "7": [0, 2, 12], "8": [0, 12], "421335633880908": 0, "h2f": 0, "4": [0, 6, 7, 8, 9, 12], "7948251330114304": 0, "f2h": 0, "30083016696826936": 0, "zeta2f": 0, "see": [0, 9, 12, 13], "battin": [0, 6, 12], "an": [0, 4, 7, 9, 10, 11, 13, 16], "introduct": 0, "mathemat": 0, "method": [0, 11, 13, 16], "astrodynam": [0, 11, 12, 15], "definit": 0, "treatment": 0, "result": [0, 10, 12, 16], "equat": [0, 7, 13, 16], "3290929552114266": 0, "f2zeta": 0, "3": [0, 2, 4, 7, 12, 13, 16], "36923933496389816": 0, "m2e_v": [0, 7], "numpi": [0, 6, 7, 8, 10, 12, 13, 14, 16], "ndarrai": [0, 13, 14, 16], "np": [0, 6, 7, 8, 10, 12, 13, 14], "linspac": [0, 7, 12], "100": [0, 7, 12], "375": 0, "shape": 0, "e2m_v": 0, "86345": 0, "m2f_v": 0, "f2m_v": [0, 7], "e2f_v": 0, "0256": 0, "f2e_v": 0, "23": 0, "n2h_v": 0, "h2n_v": 0, "n2f_v": 0, "f2n_v": 0, "h2f_v": 0, "f2h_v": 0, "zeta2f_v": 0, "f2zeta_v": 0, "design": [1, 5], "maxim": 1, "its": [1, 5, 6, 7, 9, 10, 11, 13, 15, 16], "usabl": 1, "let": [1, 7, 9, 10, 11, 12], "u": [1, 7, 9, 10, 11, 12], "know": [1, 7, 13], "what": [1, 10], "you": 1, "think": 1, "about": [1, 7, 10, 13], "global": [1, 5], "constant": [1, 5, 13], "anomali": [1, 3, 5, 12, 15, 16], "convers": [1, 5, 13], "normal": 1, "vector": [1, 6, 7, 8, 11, 12, 13, 14], "orbit": [1, 5, 7, 11, 12, 16], "element": [1, 5, 11, 13, 16], "el_typ": [1, 3, 13, 16], "ic2par": [1, 3, 13], "par2ic": [1, 3], "ic2eq": [1, 3], "eq2ic": [1, 3], "eq2par": [1, 3], "par2eq": [1, 3], "epoch": [1, 5, 10, 11, 13, 15, 16], "class": [1, 3, 5, 9, 10, 11, 14, 16], "lambert": [1, 5], "lambert_problem": [1, 6], "planet": [1, 5, 8, 10, 11], "user": [1, 2, 3, 4, 5, 9, 10, 11, 13], "implement": [1, 5, 11, 12, 13], "udpla": [1, 5, 8, 10, 11, 13], "null_udpla": [1, 16], "keplerian": [1, 3, 11, 12, 13, 16], "jpl_lp": [1, 13, 16], "tle": [1, 8, 13, 16], "spice": [1, 5, 8, 13, 15, 16], "de440": [1, 8, 10, 16], "numer": [1, 5, 10, 12], "propag": [1, 5, 11, 13, 15, 16], "propagate_lagrangian": [1, 12, 14], "access": [2, 11, 13], "number": [2, 3, 4, 6, 7, 9, 11, 12], "common": [2, 10, 11], "ar": [2, 3, 7, 9, 10, 11, 12, 13, 16], "provid": [2, 3, 5, 10, 11, 13, 15, 16], "conveni": [2, 10, 11, 12], "The": [2, 3, 4, 6, 11, 12, 13, 15, 16], "overwrit": 2, "valu": [2, 13], "need": [2, 6, 9, 10, 11], "These": [2, 3, 11, 12], "us": [2, 3, 4, 5, 7, 9, 10, 11, 12, 13, 15, 16], "intern": [2, 13], "thei": [2, 11, 13, 16], "onli": [2, 10, 11, 12], "instanti": [2, 3, 9, 10, 11, 13, 16], "object": [2, 4, 9, 10, 11, 13, 16], "": [2, 6, 7, 10, 12, 13], "unit": [2, 6, 11, 12, 13, 16], "astronom": 2, "au": [2, 11], "149597870700": 2, "cavendish": 2, "frac": [2, 3, 6, 12, 13, 14], "kg": 2, "36687e": 2, "sun": [2, 10], "gravit": [2, 6, 13, 14, 16], "paramet": [2, 3, 6, 10, 11, 13, 14, 16], "mu_sun": 2, "sec": [2, 13], "32712440018e": 2, "20": 2, "earth": [2, 10, 16], "mu_earth": 2, "398600441800000": 2, "veloc": [2, 3, 6, 10, 11, 12, 13, 14, 16], "earth_veloc": 2, "29784": 2, "691831696804": 2, "radiu": [2, 10, 11, 13, 16], "earth_radiu": 2, "6378137": 2, "j_2": 2, "earth_j2": 2, "00108262668": 2, "second": [2, 4, 6, 7, 8, 9, 11, 12, 13, 16], "dai": [2, 4, 9], "day2sec": [2, 8], "86400": 2, "degre": 2, "radian": 2, "rad2deg": 2, "57": 2, "29577951308232": 2, "default": [3, 4, 6, 9, 13, 14, 16], "oscul": [3, 13, 16], "classic": 3, "set": 3, "omega": [3, 11], "togeth": 3, "cartesian": [3, 6, 11, 12, 13, 14], "posit": [3, 6, 10, 11, 12, 13, 14, 16], "mathbf": [3, 12], "r": [3, 8, 10, 11, 12, 13, 14, 16], "v": [3, 8, 10, 11, 12, 13, 14, 16], "support": [3, 5, 9], "given": 3, "also": [3, 4, 7, 9, 12, 13], "well": [3, 5, 10, 11, 12], "equinocti": 3, "defin": [3, 4, 7, 10, 11, 12, 13, 16], "left": [3, 12], "begin": [3, 12], "arrai": [3, 6, 7, 8, 12, 13, 14, 16], "p": 3, "co": 3, "g": [3, 7, 12, 13], "sin": [3, 7], "tan": 3, "i2": 3, "right": [3, 12], "k": 3, "end": [3, 7, 11, 12], "avoid": [3, 7, 10, 13], "singular": 3, "except": [3, 13, 16], "which": [3, 7, 10, 11, 12, 13, 16], "case": [3, 7, 10, 11, 12, 13], "retrogad": 3, "convent": [3, 15], "hyperbola": 3, "enforc": 3, "thu": [3, 10, 11, 12, 15], "abl": [3, 10, 13], "where": [3, 12], "A": [3, 7, 10, 11, 16], "member": [3, 4, 6], "kep_m": 3, "kep_f": [3, 13, 16], "meq": 3, "modifi": [3, 4, 9, 13, 16], "meq_r": 3, "retrograd": [3, 6], "posvel": [3, 16], "represent": [4, 9, 11, 13], "specif": [4, 9, 12], "point": [4, 6, 9, 13], "time": [4, 6, 7, 8, 9, 10, 12, 14], "futur": [4, 9], "past": [4, 9], "rather": [4, 9, 10, 12], "confus": [4, 9, 10], "opt": [4, 9], "offer": [4, 9, 10, 11], "dedic": [4, 7, 9, 10, 12], "call": [4, 13, 16], "simpl": [4, 7, 9, 12, 13], "interfac": [4, 9, 11, 15, 16], "under": [4, 9, 10], "hood": [4, 9, 10], "seamlessli": [4, 9], "both": [4, 9, 10], "c": [4, 9, 10, 11, 13], "std": [4, 9], "chrono": [4, 9], "librari": [4, 5, 9], "python": [4, 9, 11, 13, 16], "datetim": 4, "modul": [4, 9, 10, 11, 12, 16], "julian": [4, 13, 15, 16], "date": [4, 13, 15, 16], "repres": [4, 9, 10, 13, 16], "sinc": [4, 10, 12], "start": [4, 7, 10, 11, 12, 15], "2000": [4, 9, 11, 16], "doe": [4, 9, 10, 11, 13], "account": [4, 9], "leap": [4, 9], "If": [4, 9, 10, 13], "wish": [4, 9], "exact": [4, 9], "iso": [4, 9], "8601": [4, 9], "some": [4, 9, 10, 12, 13, 15], "includ": [4, 9, 10], "he": [4, 9], "have": [4, 7, 9, 10, 12, 13], "offset": [4, 9], "himself": [4, 9], "As": [4, 9], "2023": [4, 9, 10, 16], "thi": [4, 9, 10, 11, 12, 13, 16], "mai": [4, 7, 9, 10, 13, 15], "maximum": [4, 6, 9], "28": [4, 9], "more": [4, 5, 9, 10, 11, 12, 13], "info": [4, 9, 10, 11, 13, 16], "when": [4, 8, 9, 10, 11, 13, 15], "julian_typ": [4, 9], "mjd2000": [4, 8, 9, 11, 13, 16], "construct": [4, 6, 9, 11, 13, 16], "refer": [4, 10, 11, 12, 13, 16], "jd": [4, 8, 9], "mjd": [4, 9], "12": [4, 9, 10, 11], "01": [4, 9, 10, 11, 12], "13t07": 4, "00": [4, 9, 11, 12], "000000": [4, 9, 11], "altern": [4, 16], "constructor": [4, 13, 16], "__init__": [4, 16], "str": [4, 13, 16], "string_format": 4, "string": [4, 9, 13], "format": [4, 16], "14t00": 4, "000001": 4, "year": [4, 9], "month": [4, 9], "13t00": 4, "properti": [4, 6, 13], "static": [4, 16], "now": [4, 9, 10, 12], "current": [4, 9, 10], "utc": [4, 9, 11], "coolbox": 5, "develop": 5, "european": 5, "space": [5, 7, 10, 13], "agenc": 5, "advanc": [5, 11, 12], "concept": 5, "team": 5, "Its": 5, "main": [5, 13, 16], "purpos": [5, 12, 13], "fast": 5, "prototyp": 5, "research": 5, "idea": [5, 10], "interplanetari": 5, "trajectori": [5, 11, 12], "At": [5, 7], "core": 5, "effici": [5, 12, 13], "algorithm": [5, 12], "solv": [5, 6, 7], "multipl": [5, 6, 12], "revolut": [5, 6], "problem": [5, 6, 13], "low": [5, 7, 16], "thrust": 5, "asteroid": [5, 10], "randezv": 5, "jpl": [5, 15, 16], "sgp4": [5, 8, 16], "heyoka": 5, "taylor": 5, "integr": [5, 10], "suit": 5, "ha": [5, 12], "been": 5, "dure": [5, 13], "differ": [5, 9, 10, 12], "optim": 5, "competit": 5, "gtoc": 5, "sever": 5, "paper": 5, "preliminari": 5, "mission": [5, 16], "scenario": 5, "argo": 5, "cubesat": 5, "phase": 5, "studi": 5, "titan": 5, "enceladu": 5, "tandem": 5, "analysi": 5, "hera": 5, "api": [5, 8, 11], "tutori": [5, 12], "basic": [5, 9, 12, 13], "r0": [6, 12, 14], "r1": [6, 12, 14], "tof": [6, 12, 14], "mu": [6, 10, 11, 12, 13, 14], "cw": 6, "fals": [6, 12, 13, 14], "max_rev": 6, "1d": [6, 14], "compon": [6, 11, 12, 14], "first": [6, 10, 11, 12, 13, 16], "x": 6, "y": 6, "z": 6, "xf": 6, "yf": 6, "zf": 6, "tot": [6, 14], "flight": [6, 7, 14], "bool": [6, 13, 14], "motion": [6, 10, 11, 12], "clockwis": 6, "comput": [6, 7, 10, 11, 13, 14, 15, 16], "consist": [6, 7, 11], "multirev": 6, "upon": [6, 16], "solut": 6, "store": [6, 12, 13], "data": [6, 8, 10], "lp": 6, "v0": [6, 12, 14], "1028493158958256e": 6, "16": 6, "0000000000000002": 6, "nmax": 6, "iter": 6, "attract": [6, 13, 16], "bodi": [6, 10, 11, 13, 16], "between": [6, 7, 9, 12, 13], "two": [6, 11], "v1": [6, 12, 14], "along": [6, 12], "curv": 6, "typic": 7, "mechan": 7, "indic": [7, 12], "them": [7, 10], "throughout": 7, "code": [7, 10, 11, 12, 13, 16], "document": [7, 13], "mostli": 7, "To": [7, 12], "keep": 7, "our": [7, 10, 12], "scheme": 7, "do": [7, 10, 12], "capit": 7, "letter": 7, "so": [7, 11], "transform": [7, 10], "must": [7, 10, 11, 13, 16], "singl": 7, "write": 7, "necessari": [7, 12, 16], "all": [7, 11, 12], "link": 7, "each": [7, 15], "other": [7, 11, 13], "through": [7, 10], "algebra": 7, "explicit": 7, "implicit": 7, "most": [7, 10, 11, 12, 15], "famou": 7, "kepler": 7, "here": [7, 10], "briefli": [7, 9], "showcas": [7, 10], "matplotlib": [7, 12], "pyplot": [7, 12], "plt": [7, 12], "consid": [7, 9, 10, 12, 13], "satellit": [7, 8, 10, 11, 16], "ellipt": 7, "relat": [7, 12], "rel": [7, 10], "print": [7, 8, 9, 10, 11, 12], "39017524962497735": 7, "same": [7, 10, 11, 13], "note": [7, 10], "subscript": 7, "fig": [7, 12], "figur": [7, 12], "figsiz": 7, "plot": [7, 12], "xlabel": 7, "ylabel": 7, "final": [7, 10, 12, 13, 14], "want": [7, 10], "comut": 7, "speed": [7, 12], "100000": [7, 12], "1e7": 7, "perf_count": [7, 8, 12], "rang": [7, 8, 12], "0f": 7, "4453041": 7, "satrec": 8, "export": 8, "spiceypi": [8, 16], "pyspic": 8, "line1": [8, 16], "25544u": 8, "98067a": 8, "19343": 8, "69339541": 8, "00001764": 8, "00000": [8, 16], "38792": 8, "9991": [8, 16], "line2": [8, 16], "25544": 8, "51": 8, "6439": 8, "211": 8, "2001": 8, "0007417": 8, "17": 8, "6667": 8, "85": 8, "6398": 8, "15": 8, "50103472202482": 8, "1000": 8, "20000": 8, "twoline2rv": 8, "pla": [8, 11, 13, 16], "t1": [8, 12], "process_tim": 8, "eph": [8, 13, 15, 16], "t2": [8, 12], "udplaeph": 8, "real": [8, 12], "5f": [8, 12], "42243": 8, "plaeph": 8, "19412": 8, "eph_v": [8, 13, 16], "plaephv": 8, "55561": 8, "2451544": 8, "jd_i": 8, "int": [8, 16], "jd_fr": 8, "sgp4p": 8, "40697": 8, "item": [8, 10], "b": [8, 12], "zip": [8, 12], "sgp4_arrai": 8, "sgp4pv": 8, "48900": 8, "summari": 8, "sgp4_v": 8, "util": [8, 10, 16], "load_spice_kernel": [8, 10, 16], "bsp": [8, 10, 16], "jupit": [8, 10], "barycent": [8, 10, 16], "eclipj2000": [8, 10, 16], "ssb": [8, 10, 16], "28086": 8, "37163": 8, "12372": 8, "rv": [8, 12], "_": 8, "spkezr": 8, "none": [8, 13], "spicep": [8, 12], "18197": 8, "spicev": 8, "spicepv": 8, "10087": 8, "spice_v": 8, "take": 9, "care": [9, 11], "show": [9, 12], "creat": [9, 10], "four": [9, 12], "wai": [9, 12], "pass": [9, 13], "histor": [9, 10], "directli": [9, 10], "request": [9, 14], "specifi": [9, 10], "othewis": 9, "context": 9, "arithmet": 9, "alwai": [9, 10, 11], "01t00": [9, 11], "durat": 9, "ep": [9, 13, 16], "screen": [9, 11], "explicitli": 9, "mention": 9, "type": [9, 11, 13, 16], "than": 9, "2460676": 9, "5000000": 9, "2025": [9, 10, 16], "correspond": [9, 10, 13], "30t21": 9, "35": 9, "05": [9, 12], "260406": 9, "28t00": 9, "02": [9, 10, 12], "120000": 9, "builtin": 9, "dt": [9, 12], "2033": 9, "11": [9, 10], "hour": 9, "minut": 9, "22": [9, 16], "microsecond": 9, "14532": 9, "12t12": 9, "014532": 9, "63913": 9, "51541683486": 9, "addit": 9, "subtract": 9, "timedelta": 9, "assum": [9, 12, 13, 16], "21": 9, "2353525": 9, "interpret": 9, "22t05": 9, "38": 9, "54": 9, "456000": 9, "comparison": 9, "oper": 9, "turn": 9, "handi": 9, "uniqu": 10, "move": [10, 11, 13], "comet": [10, 11], "whose": [10, 13], "fit": 10, "observ": [10, 16], "simul": 10, "encapsul": 10, "naif": [10, 16], "kernel": [10, 16], "avail": [10, 12, 15], "relev": 10, "system": [10, 16], "respect": [10, 12], "solar": [10, 16], "releas": [10, 15], "2022": 10, "accur": 10, "ones": 10, "download": 10, "binari": 10, "file": [10, 16], "contain": [10, 13, 16], "actual": [10, 12, 13], "distribut": 10, "skip": 10, "step": 10, "get": [10, 12, 15], "path": [10, 16], "de440s_kernel": 10, "kernel_fil": [10, 16], "home": 10, "runner": 10, "local": 10, "lib": 10, "python3": 10, "site": 10, "packag": [10, 16], "id": 10, "naifid": 10, "inspect_spice_kernel": 10, "readibl": 10, "naifid2nam": 10, "mercuri": [10, 16], "venu": 10, "mar": 10, "saturn": 10, "uranu": 10, "neptun": 10, "pluto": 10, "moon": 10, "nice": 10, "inspect": 10, "realiz": 10, "non": [10, 11, 12, 16], "barycentr": 10, "proce": 10, "task": 10, "jupyt": 10, "thing": 10, "usag": 10, "pre": [10, 16], "load": [10, 16], "memori": [10, 16], "done": 10, "onc": [10, 11, 13], "forget": 10, "unless": 10, "issu": 10, "unload": 10, "unload_spice_kernel": 10, "form": 10, "nveloc": 10, "722180808588": 10, "1804": 10, "157535374702": 10, "5074": 10, "16810696007": 10, "16372": 10, "2933": 10, "2858571285688": 10, "13378": 10, "581606366935": 10, "115066760074676": 10, "And": [10, 11], "python_udpla": [10, 11], "central": [10, 11, 13], "safe": [10, 11, 13, 16], "extra": [10, 11, 13, 16], "frame": [10, 13, 16], "how": [10, 11, 12, 13, 15], "mani": [10, 13], "physic": [10, 11], "particular": 10, "interpol": [10, 11], "tabl": 10, "henc": [10, 11], "ani": [10, 11, 13], "present": 10, "abov": [10, 12, 13], "valid": 10, "gener": [10, 13], "work": [10, 15], "rover": 10, "spacecraft": [10, 16], "clearli": 10, "backdraw": 10, "correct": [10, 16], "anywai": 10, "pattern": 10, "queri": 10, "automat": 10, "ship": [10, 16], "156005590351": 10, "0843": 10, "743270596831": 10, "1477": 10, "6573233296": 10, "777874": 10, "12935": 10, "993235030832": 10, "3306": 10, "5234815642566": 10, "275": 10, "73217606979927": 10, "either": [10, 13], "perform": [10, 12], "ourselv": 10, "rotat": 10, "matrix": [10, 12, 14], "would": [10, 12], "rotation_matrix": 10, "j2000": [10, 16], "inerti": 10, "orient": 10, "depend": [10, 13], "r_j2000": 10, "dot": 10, "56005590e": 10, "6": [10, 11, 12, 13, 14, 16], "84552122e": 10, "89625240e": 10, "obtain": 10, "jupiter_j2000": 10, "684552121902": 10, "1022": 10, "289625240455": 10, "7204": 10, "regardless": 11, "whether": [11, 13], "underli": [11, 13], "simpli": 11, "base": [11, 12, 13], "predict": 11, "unifi": 11, "eras": [11, 13], "hi": [11, 13], "own": [11, 13], "mandatori": [11, 13, 16], "treat": 11, "uniformli": 11, "upla": 11, "heterogen": 11, "techniqu": 11, "third": [11, 16], "parti": [11, 16], "appear": 11, "For": [11, 12, 16], "alreadi": 11, "spacecarft": [11, 16], "without": [11, 13], "describ": [11, 13, 16], "mu_central_bodi": [11, 16], "circular": 11, "dimension": [11, 12], "often": 11, "Of": 11, "cours": 11, "everyth": 11, "make": [11, 12], "sens": 11, "check": [11, 13], "anyth": [11, 13], "si": [11, 13, 16], "mix": 11, "kep3": 11, "semi": 11, "major": 11, "axi": 11, "684587122268445e": 11, "inclin": 11, "deg": 11, "big": 11, "small": 11, "anomli": 11, "ref": 11, "textual": 11, "part": 11, "report": 11, "eaul": 11, "instead": [11, 13], "origin": 11, "essenti": 11, "whatev": 11, "option": [11, 13, 16], "extra_info": 11, "expos": [11, 13], "state": [11, 12, 14], "One": 12, "coeffici": 12, "scalar": 12, "f_t": 12, "g_t": 12, "r_0": 12, "v_0": 12, "gt": 12, "analyt": 12, "express": 12, "term": 12, "univers": 12, "detail": [12, 13], "deriv": 12, "found": [12, 16], "semin": 12, "book": 12, "richard": 12, "fundament": 12, "propagate_lagrangian_v": 12, "kplerian": 12, "transit": [12, 14], "fail": [12, 13, 15, 16], "perfectli": 12, "parabol": 12, "notebook": [12, 15], "few": 12, "mpl_toolkit": 12, "mplot3d": 12, "inlin": 12, "profil": 12, "bit": 12, "cpu": 12, "re": 12, "10458": 12, "29819": 12, "period": [12, 13], "circulr": 12, "necessarili": 12, "faster": 12, "sure": 12, "t_grid": 12, "orbit1": 12, "orbit2": 12, "orbit3": 12, "pos1": 12, "pos2": 12, "pos3": 12, "just": [12, 13], "ax": 12, "project": 12, "3d": 12, "plot3d": 12, "grai": 12, "scatter3d": 12, "view_init": 12, "grid": 12, "off": 12, "ok": 12, "admittedli": 12, "underwhelm": 12, "veri": 12, "try": 12, "someth": 12, "ipython": 12, "displai": 12, "imag": 12, "filenam": 12, "sf_diagram": 12, "png": 12, "With": 12, "diagram": 12, "three": 12, "segment": 12, "ballist": 12, "arc": 12, "length": 12, "gradient": 12, "r_f": 12, "v_f": 12, "initi": [12, 14], "x_f": 12, "appli": 12, "delta": 12, "v_i": 12, "That": 12, "shall": 12, "inform": [12, 13], "bi": 12, "transfer": 12, "straightforward": 12, "after": [12, 13, 14], "convinc": 12, "m_": 12, "r_": 12, "v_": 12, "order": [12, 13], "variat": 12, "w": [12, 16], "t": [12, 13], "node": 12, "x_0": 12, "x_1": 12, "x_2": 12, "trivial": 12, "easili": 12, "seek": 12, "partial": 12, "v_1": 12, "i_v": 12, "introduc": 12, "select": [12, 16], "iv": 12, "llllll": 12, "shown": 12, "complex": [12, 13], "linear": 12, "matric": 12, "byth": 12, "condit": 12, "dv1": 12, "dv2": 12, "rv1": 12, "m1": 12, "add": 12, "rv2": 12, "m2": 12, "r2": 12, "v2": 12, "rv3": 12, "m3": 12, "readi": 12, "34626570e": 12, "9": 12, "18636660e": 12, "53491478e": 12, "03": [12, 16], "62743322e": 12, "67251795e": 12, "31214685e": 12, "01442012e": 12, "53986122e": 12, "22422462e": 12, "76624971e": 12, "74892799e": 12, "87677130e": 12, "04": 12, "20563364e": 12, "27842618e": 12, "28241661e": 12, "38362020e": 12, "94369911e": 12, "65801880e": 12, "10523246e": 12, "16229111e": 12, "86507297e": 12, "28398676e": 12, "63479977e": 12, "06363208e": 12, "91482558e": 12, "85525026e": 12, "56526220e": 12, "72411072e": 12, "06888007e": 12, "48374847e": 12, "30386029e": 12, "82337379e": 12, "43593389e": 12, "53366895e": 12, "49028595e": 12, "36023185e": 12, "diag": 12, "78070423e": 12, "08732311e": 12, "65602612e": 12, "26449980e": 12, "37560753e": 12, "10519096e": 12, "80792351e": 12, "12326674e": 12, "35231810e": 12, "17808758e": 12, "68304105e": 12, "15034626e": 12, "91248018e": 12, "27967105e": 12, "45214336e": 12, "37037657e": 12, "47952071e": 12, "03307087e": 12, "v_2": 12, "50686834e": 12, "20177166e": 12, "83470131e": 12, "20388630e": 12, "49820146e": 12, "45873151e": 12, "83271278e": 12, "45719332e": 12, "49496071e": 12, "01292574e": 12, "25392660e": 12, "80275429e": 12, "25474976e": 12, "97050307e": 12, "60311875e": 12, "79501356e": 12, "59713106e": 12, "90123105e": 12, "ephemerid": [13, 16], "possibli": 13, "etc": 13, "short": 13, "instanc": 13, "everi": 13, "least": 13, "def": 13, "self": 13, "expect": 13, "chosen": 13, "should": 13, "minim": 13, "could": 13, "fix": 13, "zero": 13, "get_mu_central_bodi": [13, 16], "get_mu_self": 13, "get_radiu": 13, "get_safe_radiu": 13, "elements_typ": 13, "get_nam": [13, 16], "get_extra_info": [13, 16], "rais": 13, "notimplementederror": 13, "unspecifi": 13, "thrown": 13, "udp": 13, "invok": 13, "deep": 13, "copi": 13, "failur": 13, "intersect": 13, "error": 13, "mismatch": 13, "signatur": 13, "otherwis": 13, "demand": 13, "over": 13, "companion": 13, "loop": 13, "behaviour": 13, "chang": 13, "who": 13, "len": 13, "extract": 13, "within": 13, "oppos": 13, "typeerror": 13, "_keplerian": 13, "my_udpla": 13, "pla2": 13, "p2": 13, "__main__": 13, "0x7ff68b63d210": 13, "0x7f8f7241c350": 13, "output": 13, "empti": 13, "averag": 13, "mainli": 13, "planetari": 13, "fly": 13, "manouvr": 13, "atmospher": 13, "circumv": 13, "radiat": 13, "environ": 13, "is_": 13, "compar": 13, "sqrt": 13, "els": [13, 16], "stm": 14, "x0": 14, "y0": 14, "z0": 14, "vx0": 14, "vy0": 14, "vz0": 14, "tupl": 14, "launch": 15, "onlin": 15, "interact": 15, "thank": 15, "infrastructur": 15, "binder": 15, "look": 15, "rocket": 15, "icon": 15, "top": 15, "page": 15, "featur": 15, "yet": 15, "latest": 15, "stabl": [15, 16], "might": 15, "execut": 15, "correctli": 15, "gist": 15, "deal": 15, "notat": 15, "de": 15, "lagrang": 15, "moot": 16, "elem": 16, "unkown": 16, "added_param": 16, "elem_typ": 16, "deafulet": 16, "keplrian": 16, "my_pla": 16, "velocti": 16, "precis": 16, "model": 16, "http": 16, "ssd": 16, "nasa": 16, "gov": 16, "approx_po": 16, "html": 16, "line": 16, "equinox": 16, "teme": 16, "33773u": 16, "97051l": 16, "23290": 16, "57931959": 16, "00002095": 16, "65841": 16, "33773": 16, "86": 16, "4068": 16, "33": 16, "1145": 16, "0009956": 16, "224": 16, "5064": 16, "135": 16, "5336": 16, "40043565770064": 16, "31": 16, "dimens": 16, "graviat": 16, "ref_fram": 16, "ob": 16, "via": 16, "readthedoc": 16, "io": 16, "en": 16, "eclipt": 16, "interest": 16, "data_archiv": 16, "ftp": 16, "www": 16, "cosmo": 16, "esa": 16, "web": 16, "rise": 16, "440": 16, "preload": 16, "body_list": 16, "quieri": 16, "befor": 16, "possibl": 16, "full": 16}, "objects": {"pykep": [[0, 0, 1, "", "e2f"], [0, 0, 1, "", "e2f_v"], [0, 0, 1, "", "e2m"], [0, 0, 1, "", "e2m_v"], [3, 1, 1, "", "el_type"], [4, 1, 1, "", "epoch"], [3, 0, 1, "", "eq2ic"], [3, 0, 1, "", "eq2par"], [0, 0, 1, "", "f2e"], [0, 0, 1, "", "f2e_v"], [0, 0, 1, "", "f2h"], [0, 0, 1, "", "f2h_v"], [0, 0, 1, "", "f2m"], [0, 0, 1, "", "f2m_v"], [0, 0, 1, "", "f2n"], [0, 0, 1, "", "f2n_v"], [0, 0, 1, "", "f2zeta"], [0, 0, 1, "", "f2zeta_v"], [0, 0, 1, "", "h2f"], [0, 0, 1, "", "h2f_v"], [0, 0, 1, "", "h2n"], [0, 0, 1, "", "h2n_v"], [3, 0, 1, "", "ic2eq"], [3, 0, 1, "", "ic2par"], [6, 1, 1, "", "lambert_problem"], [0, 0, 1, "", "m2e"], [0, 0, 1, "", "m2e_v"], [0, 0, 1, "", "m2f"], [0, 0, 1, "", "m2f_v"], [0, 0, 1, "", "n2f"], [0, 0, 1, "", "n2f_v"], [0, 0, 1, "", "n2h"], [0, 0, 1, "", "n2h_v"], [3, 0, 1, "", "par2eq"], [3, 0, 1, "", "par2ic"], [13, 1, 1, "", "planet"], [14, 1, 1, "", "propagate_lagrangian"], [0, 0, 1, "", "zeta2f"], [0, 0, 1, "", "zeta2f_v"]], "pykep.epoch": [[4, 2, 1, "", "jd"], [4, 1, 1, "", "julian_type"], [4, 2, 1, "", "mjd"], [4, 2, 1, "", "mjd2000"], [4, 3, 1, "", "now"], [4, 1, 1, "", "string_format"]], "pykep.lambert_problem": [[6, 2, 1, "", "Nmax"], [6, 2, 1, "", "iters"], [6, 2, 1, "", "mu"], [6, 2, 1, "", "r0"], [6, 2, 1, "", "r1"], [6, 2, 1, "", "tof"], [6, 2, 1, "", "v0"], [6, 2, 1, "", "v1"], [6, 2, 1, "", "x"]], "pykep.planet": [[13, 3, 1, "", "elements"], [13, 3, 1, "", "eph"], [13, 3, 1, "", "eph_v"], [13, 3, 1, "", "extract"], [13, 3, 1, "", "get_extra_info"], [13, 3, 1, "", "get_mu_central_body"], [13, 3, 1, "", "get_mu_self"], [13, 3, 1, "", "get_name"], [13, 3, 1, "", "get_radius"], [13, 3, 1, "", "get_safe_radius"], [13, 3, 1, "", "is_"], [13, 3, 1, "", "period"]], "pykep.udpla": [[16, 1, 1, "", "de440s"], [16, 1, 1, "", "jpl_lp"], [16, 1, 1, "", "keplerian"], [16, 1, 1, "", "null_udpla"], [16, 1, 1, "", "spice"], [16, 1, 1, "", "tle"]], "pykep.udpla.de440s": [[16, 3, 1, "", "body_list"], [16, 3, 1, "", "get_name"], [16, 3, 1, "", "kernel_file"]], "pykep.udpla.spice": [[16, 3, 1, "", "eph"], [16, 3, 1, "", "get_extra_info"], [16, 3, 1, "", "get_name"]], "pykep.udpla.tle": [[16, 3, 1, "", "eph"], [16, 3, 1, "", "eph_v"], [16, 3, 1, "", "get_extra_info"], [16, 3, 1, "", "get_mu_central_body"], [16, 3, 1, "", "get_name"]]}, "objtypes": {"0": "py:function", "1": "py:class", "2": "py:property", "3": "py:method"}, "objnames": {"0": ["py", "function", "Python function"], "1": ["py", "class", "Python class"], "2": ["py", "property", "Python property"], "3": ["py", "method", "Python method"]}, "titleterms": {"anomali": [0, 7], "convers": 0, "normal": 0, "vector": 0, "api": 1, "content": [1, 5], "global": 2, "constant": 2, "pykep": [2, 5, 7], "orbit": 3, "element": 3, "epoch": [4, 9], "class": [4, 6, 13], "welcom": 5, "": 5, "document": 5, "lambert": 6, "The": [7, 9, 10], "variou": 7, "julian": 9, "date": 9, "datetim": 9, "interoper": 9, "math": 9, "interfac": 10, "spice": 10, "jpl": 10, "de": 10, "eph": 10, "440": 10, "ephemerid": [10, 11], "lagrang": 12, "propag": [12, 14], "comput": 12, "overal": 12, "stm": 12, "multi": 12, "impuls": 12, "leg": 12, "planet": [13, 16], "numer": 14, "tutori": 15, "basic": 15, "list": 16, "user": 16, "implement": 16, "udpla": 16}, "envversion": {"sphinx.domains.c": 2, "sphinx.domains.changeset": 1, "sphinx.domains.citation": 1, "sphinx.domains.cpp": 8, "sphinx.domains.index": 1, "sphinx.domains.javascript": 2, "sphinx.domains.math": 2, "sphinx.domains.python": 3, "sphinx.domains.rst": 2, "sphinx.domains.std": 2, "sphinx.ext.intersphinx": 1, "sphinx": 57}, "alltitles": {"Anomalies Conversions": [[0, "anomalies-conversions"]], "Normal": [[0, "normal"]], "Vectorized": [[0, "vectorized"]], "API": [[1, "api"]], "Contents": [[1, null], [5, null]], "Global constants": [[2, "global-constants"]], "Pykep global constants": [[2, "id1"]], "Orbital Elements": [[3, "orbital-elements"]], "Epoch class": [[4, "epoch-class"]], "Welcome to pykep\u2019s documentation!": [[5, "welcome-to-pykep-s-documentation"]], "Lambert class": [[6, "lambert-class"]], "The various anomalies in pykep": [[7, "the-various-anomalies-in-pykep"]], "Epochs and Julian Dates": [[9, "epochs-and-julian-dates"]], "Julian dates": [[9, "julian-dates"]], "Datetime interoperability": [[9, "datetime-interoperability"]], "The epoch math": [[9, "the-epoch-math"]], "Interfacing to SPICE and JPL DE ephs": [[10, "interfacing-to-spice-and-jpl-de-ephs"]], "": [[10, "id1"]], "The DE 440 JPL Ephemerides": [[10, "the-de-440-jpl-ephemerides"]], "Ephemerides": [[11, "ephemerides"]], "Lagrange Propagation": [[12, "lagrange-propagation"]], "Computing the overall STM for multi-impulsive legs": [[12, "computing-the-overall-stm-for-multi-impulsive-legs"]], "Planet class": [[13, "planet-class"]], "Numerical Propagation": [[14, "numerical-propagation"]], "Tutorials": [[15, "tutorials"]], "Basic": [[15, "basic"]], "List of user implemented planets (UDPLAs)": [[16, "list-of-user-implemented-planets-udplas"]]}, "indexentries": {"e2f() (in module pykep)": [[0, "pykep.e2f"]], "e2f_v() (in module pykep)": [[0, "pykep.e2f_v"]], "e2m() (in module pykep)": [[0, "pykep.e2m"]], "e2m_v() (in module pykep)": [[0, "pykep.e2m_v"]], "f2e() (in module pykep)": [[0, "pykep.f2e"]], "f2e_v() (in module pykep)": [[0, "pykep.f2e_v"]], "f2h() (in module pykep)": [[0, "pykep.f2h"]], "f2h_v() (in module pykep)": [[0, "pykep.f2h_v"]], "f2m() (in module pykep)": [[0, "pykep.f2m"]], "f2m_v() (in module pykep)": [[0, "pykep.f2m_v"]], "f2n() (in module pykep)": [[0, "pykep.f2n"]], "f2n_v() (in module pykep)": [[0, "pykep.f2n_v"]], "f2zeta() (in module pykep)": [[0, "pykep.f2zeta"]], "f2zeta_v() (in module pykep)": [[0, "pykep.f2zeta_v"]], "h2f() (in module pykep)": [[0, "pykep.h2f"]], "h2f_v() (in module pykep)": [[0, "pykep.h2f_v"]], "h2n() (in module pykep)": [[0, "pykep.h2n"]], "h2n_v() (in module pykep)": [[0, "pykep.h2n_v"]], "m2e() (in module pykep)": [[0, "pykep.m2e"]], "m2e_v() (in module pykep)": [[0, "pykep.m2e_v"]], "m2f() (in module pykep)": [[0, "pykep.m2f"]], "m2f_v() (in module pykep)": [[0, "pykep.m2f_v"]], "n2f() (in module pykep)": [[0, "pykep.n2f"]], "n2f_v() (in module pykep)": [[0, "pykep.n2f_v"]], "n2h() (in module pykep)": [[0, "pykep.n2h"]], "n2h_v() (in module pykep)": [[0, "pykep.n2h_v"]], "zeta2f() (in module pykep)": [[0, "pykep.zeta2f"]], "zeta2f_v() (in module pykep)": [[0, "pykep.zeta2f_v"]], "el_type (class in pykep)": [[3, "pykep.el_type"]], "eq2ic() (in module pykep)": [[3, "pykep.eq2ic"]], "eq2par() (in module pykep)": [[3, "pykep.eq2par"]], "ic2eq() (in module pykep)": [[3, "pykep.ic2eq"]], "ic2par() (in module pykep)": [[3, "pykep.ic2par"]], "par2eq() (in module pykep)": [[3, "pykep.par2eq"]], "par2ic() (in module pykep)": [[3, "pykep.par2ic"]], "epoch (class in pykep)": [[4, "pykep.epoch"]], "epoch.julian_type (class in pykep)": [[4, "pykep.epoch.julian_type"]], "epoch.string_format (class in pykep)": [[4, "pykep.epoch.string_format"]], "jd (pykep.epoch property)": [[4, "pykep.epoch.jd"]], "mjd (pykep.epoch property)": [[4, "pykep.epoch.mjd"]], "mjd2000 (pykep.epoch property)": [[4, "pykep.epoch.mjd2000"]], "now() (pykep.epoch static method)": [[4, "pykep.epoch.now"]], "nmax (pykep.lambert_problem property)": [[6, "pykep.lambert_problem.Nmax"]], "iters (pykep.lambert_problem property)": [[6, "pykep.lambert_problem.iters"]], "lambert_problem (class in pykep)": [[6, "pykep.lambert_problem"]], "mu (pykep.lambert_problem property)": [[6, "pykep.lambert_problem.mu"]], "r0 (pykep.lambert_problem property)": [[6, "pykep.lambert_problem.r0"]], "r1 (pykep.lambert_problem property)": [[6, "pykep.lambert_problem.r1"]], "tof (pykep.lambert_problem property)": [[6, "pykep.lambert_problem.tof"]], "v0 (pykep.lambert_problem property)": [[6, "pykep.lambert_problem.v0"]], "v1 (pykep.lambert_problem property)": [[6, "pykep.lambert_problem.v1"]], "x (pykep.lambert_problem property)": [[6, "pykep.lambert_problem.x"]], "elements() (pykep.planet method)": [[13, "pykep.planet.elements"]], "eph() (pykep.planet method)": [[13, "pykep.planet.eph"]], "eph_v() (pykep.planet method)": [[13, "pykep.planet.eph_v"]], "extract() (pykep.planet method)": [[13, "pykep.planet.extract"]], "get_extra_info() (pykep.planet method)": [[13, "pykep.planet.get_extra_info"]], "get_mu_central_body() (pykep.planet method)": [[13, "pykep.planet.get_mu_central_body"]], "get_mu_self() (pykep.planet method)": [[13, "pykep.planet.get_mu_self"]], "get_name() (pykep.planet method)": [[13, "pykep.planet.get_name"]], "get_radius() (pykep.planet method)": [[13, "pykep.planet.get_radius"]], "get_safe_radius() (pykep.planet method)": [[13, "pykep.planet.get_safe_radius"]], "is_() (pykep.planet method)": [[13, "pykep.planet.is_"]], "period() (pykep.planet method)": [[13, "pykep.planet.period"]], "planet (class in pykep)": [[13, "pykep.planet"]], "propagate_lagrangian (class in pykep)": [[14, "pykep.propagate_lagrangian"]], "body_list() (pykep.udpla.de440s method)": [[16, "pykep.udpla.de440s.body_list"]], "de440s (class in pykep.udpla)": [[16, "pykep.udpla.de440s"]], "eph() (pykep.udpla.spice method)": [[16, "pykep.udpla.spice.eph"]], "eph() (pykep.udpla.tle method)": [[16, "pykep.udpla.tle.eph"]], "eph_v() (pykep.udpla.tle method)": [[16, "pykep.udpla.tle.eph_v"]], "get_extra_info() (pykep.udpla.spice method)": [[16, "pykep.udpla.spice.get_extra_info"]], "get_extra_info() (pykep.udpla.tle method)": [[16, "pykep.udpla.tle.get_extra_info"]], "get_mu_central_body() (pykep.udpla.tle method)": [[16, "pykep.udpla.tle.get_mu_central_body"]], "get_name() (pykep.udpla.de440s method)": [[16, "pykep.udpla.de440s.get_name"]], "get_name() (pykep.udpla.spice method)": [[16, "pykep.udpla.spice.get_name"]], "get_name() (pykep.udpla.tle method)": [[16, "pykep.udpla.tle.get_name"]], "jpl_lp (class in pykep.udpla)": [[16, "pykep.udpla.jpl_lp"]], "keplerian (class in pykep.udpla)": [[16, "pykep.udpla.keplerian"]], "kernel_file() (pykep.udpla.de440s method)": [[16, "pykep.udpla.de440s.kernel_file"]], "null_udpla (class in pykep.udpla)": [[16, "pykep.udpla.null_udpla"]], "spice (class in pykep.udpla)": [[16, "pykep.udpla.spice"]], "tle (class in pykep.udpla)": [[16, "pykep.udpla.tle"]]}})
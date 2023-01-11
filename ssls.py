# -*- coding: utf-8 -*-
"""
SSLS - Star System Lightcurve Simulator
The SSLS calculates the movements and eclipses of celestial bodies and produces a video of this.<br>
Specify mass, radius and other properties of some stars and planets in a configuration file. Then run "ssls.py <Configfilename>" to produce the video.<br>
The video shows simultanously a view of the star system from the top and from the side and the lightcurve of the system's total luminosity over time.<br>
Usually you do not need to look at or even modify the python code. Instead control the program's outcome with the config file. The meaning of all program parameters is documented in the config file.<br>
SSLS uses ffmpeg to convert the data into a video. Download ffmpeg from https://www.ffmpeg.org/download.html. Extract the zip file and add "<yourdriveandpath>\FFmpeg\bin" to Environment Variable PATH.<br>
<br>
Your questions and comments are welcome.<br>
Just open an issue on https://github.com/lichtgestalter/ssls/issues to get my attention :)<br>
"""

import configparser
import math
import matplotlib
import matplotlib.animation
import matplotlib.pyplot as plt
import numpy as np
import sys
import time

# sys.argv[1]="ssls.ini"  # If you run this script in a jupyter notebook, uncomment this line in order to provide the name of your config file.


class Parameters:
    def __init__(self, config, standard_sections):
        """Read program parameters and properties of the physical bodies from config file."""

        # [Astronomical Constants]
        g = eval(config.get("Astronomical Constants", "g"))  # For ease of use of these constants in the config file they are additionally defined here without the prefix "self.".
        au = eval(config.get("Astronomical Constants", "au"))
        r_sun = eval(config.get("Astronomical Constants", "r_sun"))
        m_sun = eval(config.get("Astronomical Constants", "m_sun"))
        l_sun = eval(config.get("Astronomical Constants", "l_sun"))
        r_jup = eval(config.get("Astronomical Constants", "r_jup"))
        m_jup = eval(config.get("Astronomical Constants", "m_jup"))
        r_earth = eval(config.get("Astronomical Constants", "r_earth"))
        m_earth = eval(config.get("Astronomical Constants", "m_earth"))
        v_earth = eval(config.get("Astronomical Constants", "v_earth"))
        self.g, self.au, self.r_sun, self.m_sun, self.l_sun, self.r_jup, self.m_jup, self.r_earth, self.m_earth, self.v_earth = g, au, r_sun, m_sun, l_sun, r_jup, m_jup, r_earth, m_earth, v_earth

        # [Video]
        self.video_file = config.get("Video", "video_file")
        self.frames = eval(config.get("Video", "frames"))
        self.fps = eval(config.get("Video", "fps"))
        self.dt = eval(config.get("Video", "dt"))
        self.sampling_rate = eval(config.get("Video", "sampling_rate"))
        self.iterations = self.frames * self.sampling_rate

        # [Scale]
        self.scope_ecl = eval(config.get("Scale", "scope_ecl"))
        self.star_scale_ecl = eval(config.get("Scale", "star_scale_ecl"))
        self.planet_scale_ecl = eval(config.get("Scale", "planet_scale_ecl"))
        self.scope_top = eval(config.get("Scale", "scope_top"))
        self.star_scale_top = eval(config.get("Scale", "star_scale_top"))
        self.planet_scale_top = eval(config.get("Scale", "planet_scale_top"))

        # [Plot]
        self.figure_width = eval(config.get("Plot", "figure_width"))
        self.figure_height = eval(config.get("Plot", "figure_height"))
        self.xlim = eval(config.get("Plot", "xlim"))
        self.ylim = eval(config.get("Plot", "ylim"))
        self.time_units = {"s": 1, "min": 60, "h": 3600, "d": 24 * 3600, "mon": 365.25 * 24 * 3600 / 12, "y": 365.25 * 24 * 3600}
        self.x_unit_name = config.get("Plot", "x_unit")
        self.x_unit_value = self.time_units[self.x_unit_name]
        self.red_dot_height = eval(config.get("Plot", "red_dot_height"))
        self.red_dot_width = eval(config.get("Plot", "red_dot_width"))

        # Checking all parameters defined so far
        for key in vars(self):
            if type(getattr(self, key)) not in [str, dict]:
                if getattr(self, key) <= 0:
                    raise Exception(f"No parameter in sections {standard_sections} may be zero or negative.")


class Body:
    def __init__(self, name, body_type, mass, radius, luminosity, startposition, velocity, a, e, i, Ω, ω, ϖ, L, ma, ea, nu, T, t, main_gravity_body, beta, color):
        """Initialize instance of physical body."""
        g, au, r_sun, m_sun, l_sun, r_jup, m_jup, r_earth, m_earth, v_earth = P.g, P.au, P.r_sun, P.m_sun, P.l_sun, P.r_jup, P.m_jup, P.r_earth, P.m_earth, P.v_earth  # For ease of use of these constants in the config file they are additionally defined here without the prefix "p.".
        self.name = name  # name
        self.body_type = body_type  # "star"/"planet"/"barycenter"
        self.mass = mass  # [kg]
        self.radius = radius  # [m]
        self.area_2d = math.pi * radius ** 2  # [m**2]
        self.luminosity = luminosity  # [W]
        self.brightness = luminosity / self.area_2d  # luminosity per (apparent) area [W/m**2]
        self.positions = np.zeros((P.iterations, 3), dtype=float)  # position for each frame
        self.main_gravity_body = main_gravity_body  # Used for calculating mu
        self.color = color  # (R, G, B)  each between 0 and 1

        if body_type == "planet":
            self.a = a  # [m] semi-major axis
            self.e = e  # [1] eccentricity
            self.i = None if i is None else math.radians(i)  # [deg] inclination
            self.Ω = None if Ω is None else math.radians(Ω)  # [deg] longitude of ascending node
            self.ω = None if ω is None else math.radians(ω)  # [deg] argument of periapsis
            self.ϖ = None if ϖ is None else math.radians(ϖ)  # [deg] longitude of periapsis
            self.L = None if L is None else math.radians(L)  # [deg] mean longitude
            self.ma = None if ma is None else math.radians(ma)  # [deg] mean anomaly
            self.ea = None if ea is None else math.radians(ea)  # [deg] eccentric anomaly
            self.nu = None if nu is None else math.radians(nu)  # [deg] true anomaly. Per definition = 270° at the time of an exoplanet's primary transit.
            self.T = T  # [s] Time of periapsis
            self.t = t  # [s] time since last time of transit
            self.ma, self.ea, self.T = None, None, None  # [rad] Only true anomaly or mean_anomaly or eccentric_anomaly or time_of_periapsis has to be provided.
            self.mu = None  # Gravitational Parameter. Depends on the masses of at least 2 bodies.
            extrascale_ecl, extrascale_top = P.planet_scale_ecl, P.planet_scale_top  # Scale radius in plot.
            self.beta = None  # unnecessary line of code?

        if body_type == "star":
            extrascale_ecl, extrascale_top = P.star_scale_ecl, P.star_scale_top  # It's a star. Scale radius in plot accordingly.
            self.beta = beta  # [1] limb darkening

        if startposition is not None and velocity is not None:  # State vectors are already in config file.
            pos = []
            for x in startposition.split(","):
                pos.append(eval(x))
            vel = []
            for x in velocity.split(","):
                vel.append(eval(x))
            self.positions[0] = np.array(pos, dtype=float)  # [m] initial position
            self.velocity = np.array(vel, dtype=float)  # [m/s]
        else:  # State vectors are not in config file. They will be calculated from Kepler orbit parameters later on after all bodies are initialized.
            self.velocity = None

        self.circle_top = matplotlib.patches.Circle((0, 0), radius * extrascale_top / P.scope_top)  # Matplotlib patch for top view
        self.circle_ecl = matplotlib.patches.Circle((0, 0), radius * extrascale_ecl / P.scope_ecl)  # Matplotlib patch for eclipsed view

        self.d, self.h, self.angle, self.eclipsed_area = 0.0, 0.0, 0.0, 0.0  # Used for calculation of eclipsed area in function eclipsed_by.

    def keplerian_elements_to_state_vectors(self):
        """Calculates the state vectors (position and velocity) from Keplerian Orbit Elements.
        Returns also true anomaly, eccentric anomaly, mean anomaly and the time of periapsis.
        [a]: https://web.archive.org/web/20160418175843/https://ccar.colorado.edu/asen5070/handouts/cart2kep2002.pdf
        [b]: https://web.archive.org/web/20170810015111/http://ccar.colorado.edu/asen5070/handouts/kep2cart_2002.doc
        [c]: https://space.stackexchange.com/questions/19322/converting-orbital-elements-to-cartesian-state-vectors/19335#19335
        [d]: https://space.stackexchange.com/questions/55356/how-to-find-eccentric-anomaly-by-mean-anomaly
        [e]: https://github.com/alfonsogonzalez/AWP/blob/main/src/python_tools/numerical_tools.py
        Numbers in comments refer to numbered formulas in [a] and [b].
        Code based on [c]. Added calculation of eccentric anomaly based on the explanations in [d] using a stripped down version of [e]."""
        a, e, i, Ω, ω, ϖ, L, ma, ea, nu, T, t, mu = self.a, self.e, self.i, self.Ω, self.ω, self.ϖ, self.L, self.ma, self.ea, self.nu, self.T, self.t, self.mu  # for readability of the following formulas

        if ω is None and ϖ is not None and Ω is not None:
            ω = ϖ - Ω
        if ma is None and L is not None and ϖ is not None:
            ma = L - ϖ
        if ea is not None:  # ea provided
            nu = 2 * math.atan(math.sqrt((1 + e) / (1 - e)) * math.tan(ea / 2))  # 3b: true anomaly (from eccentric anomaly)
            ma = ea - e * math.sin(ea)  # 2b: Mean anomaly (from eccentric anomaly). Just for completeness.
        else:  # ea not provided
            if nu is not None:  # nu provided
                ea = 2 * math.atan(math.sqrt((1 - e) / (1 + e)) * math.tan(nu / 2))  # 11a: eccentric anomaly (from true anomaly) [rad]
                ma = ea - e * math.sin(ea)  # 2b: Mean anomaly (from eccentric anomaly). Just for completeness.
            else:  # nu, ea not provided
                if ma is not None:  # ma provided
                    ea = keplers_equation_root(e, ma, ea_guess=ma)  # A good guess is important. With guess=0 the root finder very often does not converge.
                    nu = 2 * math.atan(math.sqrt((1 + e) / (1 - e)) * math.tan(ea / 2))  # 3b: true anomaly (from eccentric anomaly)
                else:  # nu, ea, ma not provided
                    if T is not None:  # T provided
                        n = math.sqrt(mu / a ** 3)  # 1b: Mean angular motion. Not needed in this function. (Except for ma, which is not needed.)
                        ma = n * T  # 1b: Mean anomaly at time of periapsis (from angular motion).
                        ea = keplers_equation_root(e, ma, ea_guess=ma)  # A good guess is important. With guess=0 the root finder very often does not converge.
                        nu = 2 * math.atan(math.sqrt((1 + e) / (1 - e)) * math.tan(ea / 2))  # 3b: true anomaly (from eccentric anomaly)
                    else:  # nu, ea, ma, T not provided
                        raise Exception("nu or ma or ea or T has to be provided to keplerian_elements_to_state_vectors()")
        n = math.sqrt(mu / a ** 3)  # 12a: mean angular motion
        T = ma / n  # Time of periapsis (from mean anomaly and angular motion). Just for completeness.

        # Now update ma, ea and nu for a delay
        # print(f'@Transit: {math.degrees(nu) =   :4.0f}   {math.degrees(ma) =   :4.0f}   {math.degrees(ea) =   :4.0f}')
        ma += t * n  # 1b
        ea = keplers_equation_root(e, ma, ea_guess=ma)  # A good guess is important. With guess=0 the root finder very often does not converge.
        nu = 2 * math.atan(math.sqrt((1 + e) / (1 - e)) * math.tan(ea / 2))  # 3b: true anomaly (from eccentric anomaly)
        # print(f' delayed: {math.degrees(nu) =   :4.0f}   {math.degrees(ma) =   :4.0f}   {math.degrees(ea) =   :4.0f}')
        # nu = nu % (2*math.pi)
        # print(f'@Transit: {math.degrees(nu) =   :4.0f}   {math.degrees(ma) =   :4.0f}   {math.degrees(ea) =   :4.0f}')
        r = a * (1 - e * math.cos(ea))  # 4b: radius r
        h = math.sqrt(mu * a * (1 - e ** 2))  # 5b: specific angular momentum h
        x = r * (math.cos(Ω) * math.cos(ω + nu) - math.sin(Ω) * math.sin(ω + nu) * math.cos(i))  # 6b: position component x
        y = r * (math.sin(Ω) * math.cos(ω + nu) + math.cos(Ω) * math.sin(ω + nu) * math.cos(i))  # 6b: position component y
        z = r * (math.sin(i) * math.sin(ω + nu))  # 6b: position component z
        p = a * (1 - e ** 2)  # 7b: Semi-latus rectum. Used in velocity calculation.
        dx = (x * h * e / (r * p)) * math.sin(nu) - (h / r) * (math.cos(Ω) * math.sin(ω + nu) + math.sin(Ω) * math.cos(ω + nu) * math.cos(i))  # 7b: velocity component x
        dy = (y * h * e / (r * p)) * math.sin(nu) - (h / r) * (math.sin(Ω) * math.sin(ω + nu) - math.cos(Ω) * math.cos(ω + nu) * math.cos(i))  # 7b: velocity component y
        dz = (z * h * e / (r * p)) * math.sin(nu) + (h / r) * (math.cos(ω + nu) * math.sin(i))  # 7b: velocity component z
        return np.array([x, y, z]), np.array([dx, dy, dz]), nu, ma, ea, T  # state vectors

    def calc_state_vectors(self, bodies):
        """Get initial position and velocity of the physical body self."""
        other_body = body_from_name(bodies, self.main_gravity_body)
        self.mu = gravitational_parameter(self, other_body)
        if self.velocity is None:  # State vectors are not in config file. So they will be calculated from Kepler orbit parameters instead.
            pos, vel, *_ = self.keplerian_elements_to_state_vectors()
            self.positions[0] = np.array(pos, dtype=float)  # [m] initial position
            self.velocity = np.array(vel, dtype=float)  # [m/s] initial velocity

    def eclipsed_by(self, body, iteration):
        """Returns area, relative_radius
        area: Area of self which is eclipsed by body.
        relative_radius: The distance of the approximated center of the eclipsed area from the center of self as a percentage of self.radius (used for limb darkening)."""
        if body.positions[iteration][1] < self.positions[iteration][1]:  # Is body nearer to viewpoint than self? (i.e. its position has a smaller y-coordinate)
            d = distance_2d_ecl(body, self, iteration)
            if d < self.radius + body.radius:  # Does body eclipse self?
                if d <= abs(self.radius - body.radius):  # Annular (i.e. ring) eclipse or total eclipse
                    if self.radius < body.radius:  # Total eclipse
                        area = self.area_2d
                        relative_radius = 0
                        # print(f'  total: {iteration:7d}  rel.area: {area/self.area_2d*100:6.0f}%  rel.r: {relative_radius*100:6.0f}%')
                        return area, relative_radius
                    else:  # Annular (i.e. ring) eclipse
                        area = body.area_2d
                        relative_radius = d / self.radius
                        # print(f'   ring: {iteration:7d}  rel.area: {area / self.area_2d * 100:6.0f}%  rel.r: {relative_radius * 100:6.0f}%')
                        return area, relative_radius
                else:  # Partial eclipse
                    # Eclipsed area is the sum of a circle segment of self plus a circle segment of body
                    # https://de.wikipedia.org/wiki/Kreissegment  https://de.wikipedia.org/wiki/Schnittpunkt#Schnittpunkte_zweier_Kreise
                    self.d = (self.radius ** 2 - body.radius ** 2 + d ** 2) / (2 * d)  # Distance of center from self to radical axis
                    body.d = (body.radius ** 2 - self.radius ** 2 + d ** 2) / (2 * d)  # Distance of center from body to radical axis
                    body.h = body.radius + self.d - d  # Height of circle segment
                    self.h = self.radius + body.d - d  # Height of circle segment
                    body.angle = 2 * math.acos(1 - body.h / body.radius)  # Angle of circle segment
                    self.angle = 2 * math.acos(1 - self.h / self.radius)  # Angle of circle segment
                    body.eclipsed_area = body.radius ** 2 * (body.angle - math.sin(body.angle)) / 2  # Area of circle segment
                    self.eclipsed_area = self.radius ** 2 * (self.angle - math.sin(self.angle)) / 2  # Area of circle segment
                    area = body.eclipsed_area + self.eclipsed_area  # Eclipsed area is sum of two circle segments.
                    relative_radius = (self.radius + self.d - body.h) / (2 * self.radius)  # Relative distance between approximated center C of eclipsed area and center of self
                    # print(f'partial: {iteration:7d}  rel.area: {area/self.area_2d*100:6.0f}%  rel.r: {relative_radius*100:6.0f}%')
                    return area, relative_radius
            else:  # No eclipse because, seen from viewer, the bodies are not close enough to each other
                return 0.0, 0.0
        else:  # body cannot eclipse self, because self is nearer to viewer than body
            return 0.0, 0.0


def find_and_check_config_file(default):
    """Check program parameters and extract config file name from them.
    Check if config file can be opened and contains all standard sections."""
    # Check program parameters and extract config file name from them.
    if len(sys.argv) == 1:
        configfilename = default
        print(f'Using default config file {configfilename}. Specify config file name as program parameter if you want to use another config file.')
    elif len(sys.argv) == 2:
        configfilename = sys.argv[1]
        print(f'Using {configfilename} as config file.')
    else:
        configfilename = sys.argv[1]
        print(f'Using {configfilename} as config file. Further program parameters are ignored.')
    config = configparser.ConfigParser(inline_comment_prefixes='#')  # Read config file.
    if len(config.read(configfilename)) < 1:  # Can the config file be opened?
        raise Exception("""config file not found. Check program parameter. If you run this script in a jupyter notebook, add sys.argv[1]='ssls.ini' to the script. """)
    for section in Standard_sections:  # Does the config file contain all standard sections?
        if section not in config.sections():
            raise Exception(f'Section {section} missing in config file.')
    return configfilename


def init_bodies(configfilename, standard_sections):
    """Read program parameters and properties of the physical bodies from config file."""
    g, au, r_sun, m_sun, l_sun, r_jup, m_jup, r_earth, m_earth, v_earth = P.g, P.au, P.r_sun, P.m_sun, P.l_sun, P.r_jup, P.m_jup, P.r_earth, P.m_earth, P.v_earth  # For ease of use of these constants in the config file they are additionally defined here without the prefix "P.".
    bodies = []
    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(configfilename)  # Read config file. Has already been done moments before for reading the program parameters.

    # Physical bodies
    for section in config.sections():
        if section not in standard_sections:  # This section must describe a physical object.
            # ω, ϖ, L, ma, ea, nu, T, t
            bodies.append(Body(name=section,
                               body_type=config.get(section, "body_type"),
                               mass=eval(config.get(section, "mass")),
                               radius=eval(config.get(section, "radius")),
                               luminosity=None if config.get(section, "luminosity", fallback=None) is None else eval(config.get(section, "luminosity")),
                               startposition=config.get(section, "startposition", fallback=None),
                               velocity=config.get(section, "velocity", fallback=None),
                               a=None if config.get(section, "a", fallback=None) is None else eval(config.get(section, "a")),  # a bit long-winded because eval() cannot handle None
                               e=None if config.get(section, "e", fallback=None) is None else eval(config.get(section, "e")),
                               i=None if config.get(section, "i", fallback=None) is None else eval(config.get(section, "i")),
                               Ω=None if config.get(section, "longitude_of_ascending_node", fallback=None) is None else eval(config.get(section, "longitude_of_ascending_node")),
                               ω=None if config.get(section, "argument_of_periapsis", fallback=None) is None else eval(config.get(section, "argument_of_periapsis")),
                               ϖ=None if config.get(section, "longitude_of_periapsis", fallback=None) is None else eval(config.get(section, "longitude_of_periapsis")),
                               L=None if config.get(section, "L", fallback=None) is None else eval(config.get(section, "L")),
                               ma=None if config.get(section, "ma", fallback=None) is None else eval(config.get(section, "ma")),
                               ea=None if config.get(section, "ea", fallback=None) is None else eval(config.get(section, "ea")),
                               nu=None if config.get(section, "nu", fallback=None) is None else eval(config.get(section, "nu")),
                               T=None if config.get(section, "T", fallback=None) is None else eval(config.get(section, "T")),
                               t=None if config.get(section, "t", fallback=None) is None else eval(config.get(section, "t")),
                               main_gravity_body=config.get(section, "main_gravity_body", fallback=None),
                               beta=eval(config.get(section, "beta")),
                               color=tuple([eval(x) for x in config.get(section, "color").split(",")])))
    # Checking parameters of physical bodies
    if len(bodies) < 1:
        raise Exception("No physical bodies specified.")
    for body in bodies:
        if body.radius <= 0:
            raise Exception(f'{body.name} has invalid radius {body.radius}.')
        if body.mass <= 0:
            raise Exception(f'{body.name} has invalid mass {body.mass}.')
        if body.luminosity < 0:
            raise Exception(f'{body.name} has invalid luminosity {body.luminosity}.')
        if body.luminosity > 0 >= body.beta:  # if body.luminosity > 0 and body.beta <= 0:
            raise Exception(f'{body.name} has invalid limb darkening parameter beta {body.beta}.')
        for c in body.color:
            if c < 0 or c > 1:
                raise Exception(f'{body.name} has invalid color value {c}.')
    return bodies


def body_from_name(bodies, bodyname):
    """Returns the body from the list bodies, which has the name bodyname."""
    if bodyname is None:
        return None
    else:
        for body in bodies:
            if body.name == bodyname:
                return body
        return None


def keplers_equation(ea, e, ma):
    """ea: eccentric anomaly [rad], e: eccentricity, ma: mean anomaly [rad]"""
    if not -2*math.pi < ea < 2*math.pi:
        raise ValueError("eccentric anomaly ea must be in radians but is outside of the range ]-2π;2π[")
    if not -2*math.pi < ma < 2*math.pi:
        raise ValueError("mean anomaly ma must be in radians but is outside of the range ]-2π;2π[")
    if not 0 <= e < 1:
        raise ValueError("eccentricity e is outside of the range [0;1[")
    return ea - e * math.sin(ea) - ma


def keplers_equation_derivative(ea, e):
    """ea: eccentric anomaly [rad], e: eccentricity"""
    return 1.0 - e * math.cos(ea)


def keplers_equation_root(e, ma, ea_guess=0.0, tolerance=1e-10, max_steps=50):
    """Calculate the root of the Kepler Equation with the Newton–Raphson method.
        e: eccentricity, ma: mean anomaly [rad], ea_guess: eccentric anomaly [rad]. ea_guess=ma is a good start."""
    for n in range(max_steps):
        delta = keplers_equation(ea_guess, e, ma) / keplers_equation_derivative(ea_guess, e)
        if abs(delta) < tolerance:
            return ea_guess - delta
        ea_guess -= delta
    raise RuntimeError('Newton\'s root solver did not converge.')


def gravitational_parameter(body1, body2):
    """Calculate the gravitational parameter of two masses
    https://en.wikipedia.org/wiki/Standard_gravitational_parameter"""
    if body1 is None or body2 is None:
        return None
    else:
        return P.g * (body1.mass + body2.mass)


def distance_2d_ecl(body1, body2, i):
    """Return distance of the centers of 2 physical bodies as seen by a viewer (projection y->0)."""
    dx = body1.positions[i][0] - body2.positions[i][0]
    dz = body1.positions[i][2] - body2.positions[i][2]
    return math.sqrt((dx ** 2 + dz ** 2))


def limbdarkening(relative_radius, beta):
    """https://en.wikipedia.org/wiki/Limb_darkening
    https://de.wikipedia.org/wiki/Photosph%C3%A4re#Mitte-Rand-Verdunkelung
    Approximates the flux of a star at a point on the star seen from a very large distance.
    The point's apparent distance from the star's center is relative_radius * radius.
    Beta depends on the wavelength. Beta=2.3 is a good compromise for the spectrum of visible light."""
    if relative_radius >= 1:
        return 1 / (1 + beta)
    return (1 + beta * math.sqrt(1 - relative_radius ** 2)) / (1 + beta)


def total_luminosity(bodies, stars, iteration):
    """"Add luminosity of all stars in the system while checking for eclipses.
    Does not yet work correctly for eclipsed eclipses (three or more bodies in line of sight at the same time)."""
    luminosity = 0.0
    for star in stars:
        luminosity += star.luminosity
        for body in bodies:
            if body != star:
                eclipsed_area, relative_radius = star.eclipsed_by(body, iteration)
                if eclipsed_area != 0:
                    luminosity -= star.brightness * eclipsed_area * limbdarkening(relative_radius, star.beta)
    return luminosity


def calc_positions_eclipses_luminosity(bodies):
    """Calculate distances, forces, accelerations, velocities of the bodies for each iteration.
    The resulting body positions and the lightcurve are stored for later use in the animation.
    Body motion calculations inspired by https://colab.research.google.com/drive/1YKjSs8_giaZVrUKDhWLnUAfebuLTC-A5."""
    stars = [body for body in bodies if body.body_type == "star"]
    lightcurve = np.zeros(P.iterations)  # Initialize lightcurve.
    lightcurve[0] = total_luminosity(bodies, stars, 0)
    for iteration in range(1, P.iterations):
        for body1 in bodies:
            force = np.array([0.0, 0.0, 0.0])
            for body2 in bodies:
                if body1 != body2:
                    # Calculate distances between bodies:
                    distance_xyz = body2.positions[iteration - 1] - body1.positions[iteration - 1]
                    distance = math.sqrt(np.dot(distance_xyz, distance_xyz))
                    force_total = P.g * body1.mass * body2.mass / distance ** 2  # Use law of gravitation to calculate force acting on body.
                    # Compute the force of attraction in each direction:
                    x, y, z = distance_xyz[0], distance_xyz[1], distance_xyz[2]
                    polar_angle = math.acos(z / distance)
                    azimuth_angle = math.atan2(y, x)
                    force[0] += math.sin(polar_angle) * math.cos(azimuth_angle) * force_total
                    force[1] += math.sin(polar_angle) * math.sin(azimuth_angle) * force_total
                    force[2] += math.cos(polar_angle) * force_total
            acceleration = force / body1.mass  # Compute the acceleration in each direction.
            body1.velocity += acceleration * P.dt  # Compute the velocity in each direction.
            # Update positions:
            movement = body1.velocity * P.dt - 0.5 * acceleration * P.dt ** 2
            body1.positions[iteration] = body1.positions[iteration - 1] + movement
        lightcurve[iteration] = total_luminosity(bodies, stars, iteration)  # Update lightcurve.
        if iteration % int(round(P.iterations / 10)) == 0:  # Inform user about program's progress.
            print(f'{round(iteration / P.iterations * 10) * 10:3d}% ', end="")
    return lightcurve, bodies


def calc_physics(bodies):
    """Calculate body positions and the resulting lightcurve."""
    print(f'Producing {P.frames / P.fps:.0f} seconds long video, covering {P.dt * P.iterations / 60 / 60 / 24:5.2f} earth days. ({P.dt * P.sampling_rate * P.fps / 60 / 60 / 24:.2f} earth days per video second.)')
    print(f'Calculating {P.iterations:6d} iterations: ', end="")
    tic = time.perf_counter()
    lightcurve, bodies = calc_positions_eclipses_luminosity(bodies)
    lightcurve /= lightcurve.max(initial=None)  # Normalize flux.
    toc = time.perf_counter()
    print(f' {toc - tic:7.2f} seconds  ({P.iterations / (toc - tic):.0f} iterations/second)')
    return lightcurve, bodies


def tic_delta(scope):
    """Returns a distance between two tics on an axis so that the total number of tics on that axis is between 5 and 10."""
    if scope <= 0:  # no or constant values
        return 1
    delta = 10 ** np.floor(math.log10(scope))
    if scope / delta < 5:
        if scope / delta < 2:
            return delta / 5
        else:
            return delta / 2
    else:
        return delta


def init_plot(sampled_lightcurve):
    """Initialize the matplotlib figure containing 3 axis:
    Eclipse view (top left): projection (x,y,z) -> (x,z), order = -y.
    Top view (top right): projection (x,y,z) -> (x,y), order = z.
    Lightcurve (bottom)"""
    fig = plt.figure()
    fig.set_figwidth(P.figure_width)
    fig.set_figheight(P.figure_height)
    fig.set_facecolor("black")  # background color outside of ax_eclipse and ax_lightcurve
    buffer = 0
    fig.subplots_adjust(left=buffer, right=1.0 - buffer, bottom=buffer, top=1 - buffer)  # Positions of the subplots edges, as a fraction of the figure width.

    ax_eclipse = plt.subplot2grid(shape=(5, 2), loc=(0, 0), rowspan=4, colspan=1)
    ax_eclipse.set_xlim(-P.xlim, P.xlim)
    ax_eclipse.set_ylim(-P.ylim, P.ylim)
    ax_eclipse.set_aspect('equal')
    ax_eclipse.set_facecolor("black")  # background color
    # ax_eclipse.get_xaxis().set_visible(False)
    # ax_eclipse.get_yaxis().set_visible(False)

    ax_top = plt.subplot2grid(shape=(5, 2), loc=(0, 1), rowspan=4, colspan=1)
    ax_top.set_xlim(-P.xlim, P.xlim)
    ax_top.set_ylim(-P.ylim, P.ylim)
    ax_top.set_aspect('equal')
    ax_top.set_facecolor("black")  # background color

    ax_lightcurve = plt.subplot2grid(shape=(5, 1), loc=(4, 0), rowspan=1, colspan=1)
    ax_lightcurve.set_facecolor("black")  # background color

    ax_lightcurve.tick_params(axis='x', colors='grey')
    xmax = P.iterations * P.dt / P.x_unit_value
    ax_lightcurve.set_xlim(0, xmax)
    xvalues = [x * tic_delta(xmax) for x in range(round(xmax / tic_delta(xmax)))]
    xlabels = [f'{round(x, 4)} {P.x_unit_name}' for x in xvalues]
    ax_lightcurve.set_xticks(xvalues, labels=xlabels)

    ax_lightcurve.tick_params(axis='y', colors='grey')
    minl = Lightcurve.min(initial=None)
    maxl = Lightcurve.max(initial=None)
    if minl == maxl:
        minl *= 0.99
    scope = maxl - minl
    buffer = 0.05 * scope
    ax_lightcurve.set_ylim(minl - buffer, maxl + buffer)

    ticdelta = tic_delta(maxl - minl)
    yvalues = [1 - y * ticdelta for y in range(round(float((maxl - minl) / ticdelta)))]
    ylabels = [f'{round(100 * y, 10)} %' for y in yvalues]
    ax_lightcurve.set_yticks(yvalues, labels=ylabels)

    time_axis = np.arange(0, round(P.iterations * P.dt), round(P.sampling_rate * P.dt), dtype=float)
    time_axis /= P.x_unit_value
    ax_lightcurve.plot(time_axis, sampled_lightcurve[0:len(time_axis)], color="white")

    red_dot = matplotlib.patches.Ellipse((0, 0), P.iterations * P.dt * P.red_dot_width / P.x_unit_value, scope * P.red_dot_height)  # matplotlib patch
    red_dot.set(zorder=2)  # Dot in front of lightcurve.
    red_dot.set_color((1, 0, 0))  # red
    ax_lightcurve.add_patch(red_dot)
    plt.tight_layout()  # Automatically adjust padding horizontally as well as vertically.
    return fig, ax_top, ax_eclipse, ax_lightcurve, red_dot


def prepare_animation(bodies):
    """Initialize all matplotlib objects."""
    sampled_lightcurve = np.take(Lightcurve, range(0, P.iterations, P.sampling_rate))  # Use only some of the calculated positions for the animation because it is so slow.
    fig, ax_top, ax_eclipse, ax_lightcurve, red_dot = init_plot(sampled_lightcurve)  # Adjust constants in section [Plot] of config file to fit your screen.
    for body in bodies:  # Circles represent the bodies in the animation. Set their colors and add them to the matplotlib axis.
        body.circle_top.set_color(body.color)
        body.circle_ecl.set_color(body.color)
        ax_top.add_patch(body.circle_top)
        ax_eclipse.add_patch(body.circle_ecl)
    return fig, bodies, red_dot


def next_animation_frame(frame, bodies, red_dot):
    """Update patches. Send new circle positions to animation function.
    First parameter comes from iterator frames (a parameter of FuncAnimation).
    The other parameters are given to this function via the parameter fargs of FuncAnimation."""
    for body in bodies:  # Top view: projection (x,y,z) -> (x,y), order = z
        body.circle_top.set(zorder=body.positions[frame * P.sampling_rate][2])
        body.circle_top.center = body.positions[frame * P.sampling_rate][0] / P.scope_top, body.positions[frame * P.sampling_rate][1] / P.scope_top
    for body in bodies:  # Eclipse view: projection (x,y,z) -> (x,z), order = -y
        body.circle_ecl.set(zorder=-body.positions[frame * P.sampling_rate][1])
        body.circle_ecl.center = body.positions[frame * P.sampling_rate][0] / P.scope_ecl, body.positions[frame * P.sampling_rate][2] / P.scope_ecl
    red_dot.center = P.dt * P.sampling_rate * frame / P.x_unit_value, Lightcurve[frame * P.sampling_rate]
    if frame > 0 and frame % int(round(P.frames / 10)) == 0:  # Inform user about program's progress.
        print(f'{round(frame / P.frames * 10) * 10:3d}% ', end="")


def render_animation(bodies, red_dot):
    """Calls next_animation_frame() for each frame and saves the video."""
    print(f'Animating {P.frames:8d} frames:     ', end="")
    tic = time.perf_counter()
    anim = matplotlib.animation.FuncAnimation(Fig, next_animation_frame, fargs=(bodies, red_dot,), interval=1000 / P.fps, frames=P.frames, blit=False)
    anim.save(P.video_file, fps=P.fps, metadata={"title": " "}, extra_args=['-vcodec', 'libx264'])  # https://www.ffmpeg.org/libavcodec.html
    toc = time.perf_counter()
    print(f' {toc - tic:7.2f} seconds  ({P.frames / (toc - tic):.0f} frames/second)')
    print(f'{P.video_file} saved.')
    return bodies, red_dot


if __name__ == '__main__':
    Standard_sections = ["Astronomical Constants", "Video", "Plot", "Scale"]
    Config = configparser.ConfigParser(inline_comment_prefixes='#')
    Configfilename = find_and_check_config_file(default="ssls.ini")
    Config.read(Configfilename)
    P = Parameters(Config, Standard_sections)  # Read program parameters from config file.
    Bodies = init_bodies(Configfilename, Standard_sections)  # Read the properties of the physical bodies from the config file and write them into <bodies>, a list of all physical objects of the simulation.
    for B in Bodies:
        B.calc_state_vectors(Bodies)
    Lightcurve, Bodies = calc_physics(Bodies)  # Calculate body positions and the resulting lightcurve.
    Fig, Bodies, Red_dot = prepare_animation(Bodies)
    Bodies, Red_dot = render_animation(Bodies, Red_dot)

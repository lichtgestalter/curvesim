# -*- coding: utf-8 -*-
# SSLS - Star System Lightcurve Simulator
# The SSLS calculates the movements and eclipses of celestial bodies and produces a video of this.
# Specify mass, radius and other properties of some stars and planets in a configuration file.
# Then run "ssls.py <Configfilename>" to produce the video.
# The video shows simultanously a view of the star system from the top and from the side and
# the lightcurve of the system's total luminosity over time.
# Usually you do not need to look at or even modify the python code. Instead control the program's
# outcome with the config file. The meaning of all program parameters is documented in the config file.
# SSLS uses ffmpeg to convert the data into a video. Download ffmpeg from https://www.ffmpeg.org/download.html.
# Extract the zip file and add "<yourdriveandpath>\FFmpeg\bin" to Environment Variable PATH.<br>
#
# Your questions and comments are welcome.
# Just open an issue on https://github.com/lichtgestalter/ssls/issues to get my attention :)

import configparser
import math
import matplotlib
import matplotlib.animation
import numpy as np
import time

from cs_animation import CurveSimAnimation
from cs_parameters import CurveSimParameters, Standard_sections
from cs_physics import CurveSimPhysics


# If you run this script in a jupyter notebook, uncomment this line in order to provide the name of your config file.
# sys.argv[1]="ssls.ini"


class CurveSimBody:

    # noinspection NonAsciiCharacters,PyPep8Naming,PyUnusedLocal
    def __init__(self, name, body_type, mass, radius, luminosity, startposition, velocity, a, e, i, Ω, ω, ϖ, L, ma, ea,
                 nu, T, t, beta, color):
        """Initialize instance of physical body."""
        # For ease of use of constants in the config file they are additionally defined here without the prefix "p.".
        g, au, r_sun, m_sun, l_sun = P.g, P.au, P.r_sun, P.m_sun, P.l_sun
        r_jup, m_jup, r_earth, m_earth, v_earth = P.r_jup, P.m_jup, P.r_earth, P.m_earth, P.v_earth
        self.name = name  # name
        self.body_type = body_type  # "star" or "planet"
        self.mass = mass  # [kg]
        self.radius = radius  # [m]
        self.area_2d = math.pi * radius ** 2  # [m**2]
        self.luminosity = luminosity  # [W]
        self.brightness = luminosity / self.area_2d  # luminosity per (apparent) area [W/m**2]
        self.positions = np.zeros((P.iterations, 3), dtype=float)  # position for each frame
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
            self.beta = None  # unnecessary line of code?

        if body_type == "star":
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

        # Used for calculation of eclipsed area in function eclipsed_by.
        self.d, self.h, self.angle, self.eclipsed_area = 0.0, 0.0, 0.0, 0.0

    # noinspection NonAsciiCharacters,PyPep8Naming,PyUnusedLocal
    def keplerian_elements_to_state_vectors(self):
        """Calculates the state vectors (position and velocity) from Keplerian Orbit Elements.
        Returns also true anomaly, eccentric anomaly, mean anomaly and the time of periapsis.
        [a]: https://web.archive.org/web/20160418175843/https://ccar.colorado.edu/asen5070/handouts/cart2kep2002.pdf
        [b]: https://web.archive.org/web/20170810015111/http://ccar.colorado.edu/asen5070/handouts/kep2cart_2002.doc
        [c]: https://space.stackexchange.com/questions/19322/converting-orbital-elements-to-cartesian-state-vectors/19335#19335
        [d]: https://space.stackexchange.com/questions/55356/how-to-find-eccentric-anomaly-by-mean-anomaly
        [e]: https://github.com/alfonsogonzalez/AWP/blob/main/src/python_tools/numerical_tools.py
        Numbers in comments refer to numbered formulas in [a] and [b].
        Code based on [c]. Added calculation of eccentric anomaly based on the explanations
        in [d] using a stripped down version of [e]."""
        a, e, i, Ω, ω, ϖ, L = self.a, self.e, self.i, self.Ω, self.ω, self.ϖ, self.L  # for readability of formulas
        ma, ea, nu, T, t, mu = self.ma, self.ea, self.nu, self.T, self.t, self.mu  # for readability of formulas

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
                    ea = CurveSimPhysics.kepler_equation_root(e, ma, ea_guess=ma)  # A good guess is important. With guess=0 the root finder very often does not converge.
                    nu = 2 * math.atan(math.sqrt((1 + e) / (1 - e)) * math.tan(ea / 2))  # 3b: true anomaly (from eccentric anomaly)
                else:  # nu, ea, ma not provided
                    if T is not None:  # T provided
                        n = math.sqrt(mu / a ** 3)  # 1b: Mean angular motion. Not needed in this function. (Except for ma, which is not needed.)
                        ma = n * T  # 1b: Mean anomaly at time of periapsis (from angular motion).
                        ea = CurveSimPhysics.kepler_equation_root(e, ma, ea_guess=ma)  # A good guess is important. With guess=0 the root finder very often does not converge.
                        nu = 2 * math.atan(math.sqrt((1 + e) / (1 - e)) * math.tan(ea / 2))  # 3b: true anomaly (from eccentric anomaly)
                    else:  # nu, ea, ma, T not provided
                        raise Exception("nu or ma or ea or T has to be provided to keplerian_elements_to_state_vectors()")
        n = math.sqrt(mu / a ** 3)  # 12a: mean angular motion
        T = ma / n  # Time of periapsis (from mean anomaly and angular motion). Just for completeness.

        # Now update ma, ea and nu for a delay
        # print(f'@Transit: {math.degrees(nu) =   :4.0f}   {math.degrees(ma) =   :4.0f}   {math.degrees(ea) =   :4.0f}')
        ma += t * n  # 1b
        ea = CurveSimPhysics.kepler_equation_root(e, ma, ea_guess=ma)  # A good guess is important. With guess=0 the root finder very often does not converge.
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
        self.mu = CurveSimPhysics.gravitational_parameter(bodies, P.g)  # is the same for all bodies in the system, because they are orbiting a common barycenter
        if self.velocity is None:  # State vectors are not in config file. So they will be calculated from Kepler orbit parameters instead.
            pos, vel, *_ = self.keplerian_elements_to_state_vectors()
            self.positions[0] = np.array(pos, dtype=float)  # [m] initial position
            self.velocity = np.array(vel, dtype=float)  # [m/s] initial velocity

    def eclipsed_by(self, body, iteration):
        """Returns area, relative_radius
        area: Area of self which is eclipsed by body.
        relative_radius: The distance of the approximated center of the eclipsed area from the center of self as a percentage of self.radius (used for limb darkening)."""
        if body.positions[iteration][1] < self.positions[iteration][1]:  # Is body nearer to viewpoint than self? (i.e. its position has a smaller y-coordinate)
            d = CurveSimPhysics.distance_2d_ecl(body, self, iteration)
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


class CurveSimBodies(list):

    # noinspection PyUnusedLocal
    def __init__(self, configfilename, standard_sections):
        """Initialize instances of physical bodies.
        Read program parameters and properties of the bodies from config file.
        Initialize the circles in the animation (matplotlib patches)"""
        # For ease of use of these constants in the config file they are additionally defined here without the prefix "P.".
        g, au, r_sun, m_sun, l_sun = P.g, P.au, P.r_sun, P.m_sun, P.l_sun
        r_jup, m_jup, r_earth, m_earth, v_earth = P.r_jup, P.m_jup, P.r_earth, P.m_earth, P.v_earth
        config = configparser.ConfigParser(inline_comment_prefixes='#')
        config.read(configfilename)  # Read config file. (Has been done moments before for reading the program parameters.)
        super().__init__()

        # Physical bodies
        for section in config.sections():
            if section not in standard_sections:  # This section must describe a physical object.
                self.append(CurveSimBody(name=section,
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
                                         beta=eval(config.get(section, "beta")),
                                         color=tuple([eval(x) for x in config.get(section, "color").split(",")])))
        # Checking parameters of physical bodies
        if len(self) < 1:
            raise Exception("No physical bodies specified.")
        for body in self:
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
            body.calc_state_vectors(self)
        self.generate_patches()

    def total_luminosity(self, stars, iteration):
        """"Add luminosity of all stars in the system while checking for eclipses.
        Does not yet work correctly for eclipsed eclipses (three or more bodies in line of sight at the same time)."""
        luminosity = 0.0
        for star in stars:
            luminosity += star.luminosity
            for body in self:
                if body != star:
                    eclipsed_area, relative_radius = star.eclipsed_by(body, iteration)
                    if eclipsed_area != 0:
                        luminosity -= star.brightness * eclipsed_area * CurveSimPhysics.limbdarkening(relative_radius, star.beta)
        return luminosity

    def calc_positions_eclipses_luminosity(self):
        """Calculate distances, forces, accelerations, velocities of the bodies for each iteration.
        The resulting body positions and the lightcurve are stored for later use in the animation.
        Body motion calculations inspired by https://colab.research.google.com/drive/1YKjSs8_giaZVrUKDhWLnUAfebuLTC-A5."""
        stars = [body for body in self if body.body_type == "star"]
        lightcurve = np.zeros(P.iterations)  # Initialize lightcurve.
        lightcurve[0] = self.total_luminosity(stars, 0)
        for iteration in range(1, P.iterations):
            for body1 in self:
                force = np.array([0.0, 0.0, 0.0])
                for body2 in self:
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
            lightcurve[iteration] = self.total_luminosity(stars, iteration)  # Update lightcurve.
            if iteration % int(round(P.iterations / 10)) == 0:  # Inform user about program's progress.
                print(f'{round(iteration / P.iterations * 10) * 10:3d}% ', end="")
        return lightcurve, self

    def calc_physics(self):
        """Calculate body positions and the resulting lightcurve."""
        print(f'Producing {P.frames / P.fps:.0f} seconds long video, covering {P.dt * P.iterations / 60 / 60 / 24:5.2f} '
              f'earth days. ({P.dt * P.sampling_rate * P.fps / 60 / 60 / 24:.2f} earth days per video second.)')
        print(f'Calculating {P.iterations:6d} iterations: ', end="")
        tic = time.perf_counter()
        lightcurve, bodies = self.calc_positions_eclipses_luminosity()
        lightcurve /= lightcurve.max(initial=None)  # Normalize flux.
        toc = time.perf_counter()
        print(f' {toc - tic:7.2f} seconds  ({P.iterations / (toc - tic):.0f} iterations/second)')
        return lightcurve

    def calc_patch_radii(self):
        """If autoscaling is on, this function calculates the radii of the circles (matplotlib patches) of the animation."""
        radius_list = [body.radius for body in self]  # radii of all bodies
        # print(f'{rlist=}')
        log_list = [math.log10(i) for i in radius_list]  # log10 of all radii
        # print(f'{log_list=}')
        log_scaled_list = [i * P.min_radius / min(log_list) for i in log_list]  # scaled log10 lineary, so the smallest circle has the desired radius
        # print(f'{min_ok=}')
        exp_numerator = math.log10((P.max_radius - max(log_scaled_list)) / max(log_scaled_list))
        exp_denominator = math.log10(max(log_list) - min(log_list))
        # print(x, y, exponent)
        radius_list = [i * (1 + (j - min(log_list)) ** (exp_numerator / exp_denominator)) for i, j in zip(log_scaled_list, log_list)]  # scaled log10 exponentially, so all circles have the desired radius
        # print(f'{fertig=}')
        for body, radius in zip(self, radius_list):
            body.patch_radius = radius

    def generate_patches(self):
        """Generates the circles (matplotlib patches) of the animation."""
        if P.autoscaling:
            print("autoscaling on")
            self.calc_patch_radii()
            for body in self:
                body.circle_top = matplotlib.patches.Circle(xy=(0, 0), radius=body.patch_radius)  # Matplotlib patch for top view
                body.circle_ecl = matplotlib.patches.Circle(xy=(0, 0), radius=body.patch_radius)  # Matplotlib patch for eclipsed view
        else:
            print("autoscaling off")
            for body in self:
                if body.body_type == "planet":
                    extrascale_ecl, extrascale_top = P.planet_scale_ecl, P.planet_scale_top  # Scale radius in plot.
                else:
                    extrascale_ecl, extrascale_top = P.star_scale_ecl, P.star_scale_top  # It's a star. Scale radius in plot accordingly.
                body.circle_top = matplotlib.patches.Circle((0, 0), radius=body.radius * extrascale_top / P.scope_top)  # Matplotlib patch for top view
                body.circle_ecl = matplotlib.patches.Circle((0, 0), radius=body.radius * extrascale_ecl / P.scope_ecl)  # Matplotlib patch for eclipsed view


if __name__ == '__main__':
    Config = configparser.ConfigParser(inline_comment_prefixes='#')
    Configfilename = CurveSimParameters.find_and_check_config_file(default="ssls.ini")
    Config.read(Configfilename)
    P = CurveSimParameters(Config, Standard_sections)  # Read program parameters from config file.
    Bodies = CurveSimBodies(Configfilename, Standard_sections)  # Read the properties of the physical bodies from the config file and write them into <Bodies>, a list of all physical objects of the simulation.
    Lightcurve = Bodies.calc_physics()  # Calculate body positions and the resulting lightcurve.
    animation = CurveSimAnimation(P, Bodies, Lightcurve)

import configparser
import math
import time
import matplotlib
import matplotlib.animation
import numpy as np

from cs_body import CurveSimBody
from cs_physics import CurveSimPhysics

class CurveSimBodies(list):

    # noinspection PyUnusedLocal
    def __init__(self, p, configfilename, standard_sections):
        """Initialize instances of physical bodies.
        Read program parameters and properties of the bodies from config file.
        Initialize the circles in the animation (matplotlib patches)"""
        # For ease of use of these constants in the config file they are additionally defined here without the prefix "p.".
        g, au, r_sun, m_sun, l_sun = p.g, p.au, p.r_sun, p.m_sun, p.l_sun
        r_jup, m_jup, r_earth, m_earth, v_earth = p.r_jup, p.m_jup, p.r_earth, p.m_earth, p.v_earth
        config = configparser.ConfigParser(inline_comment_prefixes='#')
        config.read(configfilename)  # Read config file. (Has been done moments before for reading the program parameters.)
        super().__init__()  # create object by calling the constructor of class list

        # Physical bodies
        for section in config.sections():
            if section not in standard_sections:  # This section must describe a physical object.
                self.append(CurveSimBody(p=p,
                                         name=section,
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
            body.calc_state_vectors(p, self)
        self.generate_patches(p)

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

    def calc_positions_eclipses_luminosity(self, p):
        """Calculate distances, forces, accelerations, velocities of the bodies for each iteration.
        The resulting body positions and the lightcurve are stored for later use in the animation.
        Body motion calculations inspired by https://colab.research.google.com/drive/1YKjSs8_giaZVrUKDhWLnUAfebuLTC-A5."""
        stars = [body for body in self if body.body_type == "star"]
        lightcurve = np.zeros(p.iterations)  # Initialize lightcurve.
        lightcurve[0] = self.total_luminosity(stars, 0)
        for iteration in range(1, p.iterations):
            for body1 in self:
                force = np.array([0.0, 0.0, 0.0])
                for body2 in self:
                    if body1 != body2:
                        # Calculate distances between bodies:
                        distance_xyz = body2.positions[iteration - 1] - body1.positions[iteration - 1]
                        distance = math.sqrt(np.dot(distance_xyz, distance_xyz))
                        force_total = p.g * body1.mass * body2.mass / distance ** 2  # Use law of gravitation to calculate force acting on body.
                        # Compute the force of attraction in each direction:
                        x, y, z = distance_xyz[0], distance_xyz[1], distance_xyz[2]
                        polar_angle = math.acos(z / distance)
                        azimuth_angle = math.atan2(y, x)
                        force[0] += math.sin(polar_angle) * math.cos(azimuth_angle) * force_total
                        force[1] += math.sin(polar_angle) * math.sin(azimuth_angle) * force_total
                        force[2] += math.cos(polar_angle) * force_total
                acceleration = force / body1.mass  # Compute the acceleration in each direction.
                body1.velocity += acceleration * p.dt  # Compute the velocity in each direction.
                # Update positions:
                movement = body1.velocity * p.dt - 0.5 * acceleration * p.dt ** 2
                body1.positions[iteration] = body1.positions[iteration - 1] + movement
            lightcurve[iteration] = self.total_luminosity(stars, iteration)  # Update lightcurve.
            if iteration % int(round(p.iterations / 10)) == 0:  # Inform user about program's progress.
                print(f'{round(iteration / p.iterations * 10) * 10:3d}% ', end="")
        return lightcurve, self

    def calc_physics(self, p):
        """Calculate body positions and the resulting lightcurve."""
        print(f'Producing {p.frames / p.fps:.0f} seconds long video, covering {p.dt * p.iterations / 60 / 60 / 24:5.2f} '
              f'earth days. ({p.dt * p.sampling_rate * p.fps / 60 / 60 / 24:.2f} earth days per video second.)')
        print(f'Calculating {p.iterations:6d} iterations: ', end="")
        tic = time.perf_counter()
        lightcurve, bodies = self.calc_positions_eclipses_luminosity(p)
        lightcurve /= lightcurve.max(initial=None)  # Normalize flux.
        toc = time.perf_counter()
        print(f' {toc - tic:7.2f} seconds  ({p.iterations / (toc - tic):.0f} iterations/second)')
        return lightcurve

    def calc_patch_radii(self, p):
        """If autoscaling is on, this function calculates the radii of the circles (matplotlib patches) of the animation."""
        radius_list = [body.radius for body in self]  # radii of all bodies
        # print(f'{rlist=}')
        log_list = [math.log10(i) for i in radius_list]  # log10 of all radii
        # print(f'{log_list=}')
        log_scaled_list = [i * p.min_radius / min(log_list) for i in log_list]  # scaled log10 lineary, so the smallest circle has the desired radius
        # print(f'{min_ok=}')
        exp_numerator = math.log10((p.max_radius - max(log_scaled_list)) / max(log_scaled_list))
        exp_denominator = math.log10(max(log_list) - min(log_list))
        # print(x, y, exponent)
        radius_list = [i * (1 + (j - min(log_list)) ** (exp_numerator / exp_denominator)) for i, j in zip(log_scaled_list, log_list)]  # scaled log10 exponentially, so all circles have the desired radius
        # print(f'{fertig=}')
        for body, radius in zip(self, radius_list):
            body.patch_radius = radius

    def generate_patches(self, p):
        """Generates the circles (matplotlib patches) of the animation."""
        if p.autoscaling:
            print("autoscaling on")
            self.calc_patch_radii(p)
            for body in self:
                body.circle_top = matplotlib.patches.Circle(xy=(0, 0), radius=body.patch_radius)  # Matplotlib patch for top view
                body.circle_ecl = matplotlib.patches.Circle(xy=(0, 0), radius=body.patch_radius)  # Matplotlib patch for eclipsed view
        else:
            print("autoscaling off")
            for body in self:
                if body.body_type == "planet":
                    extrascale_ecl, extrascale_top = p.planet_scale_ecl, p.planet_scale_top  # Scale radius in plot.
                else:
                    extrascale_ecl, extrascale_top = p.star_scale_ecl, p.star_scale_top  # It's a star. Scale radius in plot accordingly.
                body.circle_top = matplotlib.patches.Circle((0, 0), radius=body.radius * extrascale_top / p.scope_top)  # Matplotlib patch for top view
                body.circle_ecl = matplotlib.patches.Circle((0, 0), radius=body.radius * extrascale_ecl / p.scope_ecl)  # Matplotlib patch for eclipsed view

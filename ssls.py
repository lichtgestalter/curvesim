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
import numpy as np
import time

from cs_animation import CurveSimAnimation
from cs_body import CurveSimBody
from cs_parameters import CurveSimParameters, Standard_sections
from cs_physics import CurveSimPhysics


# If you run this script in a jupyter notebook, uncomment this line in order to provide the name of your config file.
# sys.argv[1]="ssls.ini"


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
                self.append(CurveSimBody(p=P,
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
            body.calc_state_vectors(P, self)
        CurveSimBody.generate_patches(P, self)

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


if __name__ == '__main__':
    Config = configparser.ConfigParser(inline_comment_prefixes='#')
    Configfilename = CurveSimParameters.find_and_check_config_file(default="ssls.ini")
    Config.read(Configfilename)
    P = CurveSimParameters(Config, Standard_sections)  # Read program parameters from config file.
    Bodies = CurveSimBodies(Configfilename, Standard_sections)  # Read the properties of the physical bodies from the config file and write them into <Bodies>, a list of all physical objects of the simulation.
    Lightcurve = Bodies.calc_physics()  # Calculate body positions and the resulting lightcurve.
    animation = CurveSimAnimation(P, Bodies, Lightcurve)

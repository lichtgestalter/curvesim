# -*- coding: utf-8 -*-
# SSLS - Star System Lightcurve Simulator
# The SSLS calculates the movements and eclipses of celestial bodies and produces a video of this.
# Specify mass, radius and other properties of some stars and planets in a configuration file.
# Then run "cs_main.py <Configfilename>" to produce the video.
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

from cs_animation import CurveSimAnimation
from cs_bodies import CurveSimBodies
from cs_parameters import CurveSimParameters, Standard_sections

# import configparser
# import math
# import matplotlib
# import matplotlib.animation
# import numpy as np
# import time
#
# from cs_animation import CurveSimAnimation
# from cs_body import CurveSimBody
# from cs_parameters import CurveSimParameters, Standard_sections
# from cs_physics import CurveSimPhysics


# If you run this script in a jupyter notebook, uncomment this line in order to provide the name of your config file.
# sys.argv[1]="ssls.ini"


if __name__ == '__main__':
    Config = configparser.ConfigParser(inline_comment_prefixes='#')
    Configfilename = CurveSimParameters.find_and_check_config_file(default="ssls.ini")
    Config.read(Configfilename)
    P = CurveSimParameters(Config, Standard_sections)  # Read program parameters from config file.
    Bodies = CurveSimBodies(P, Configfilename, Standard_sections)  # Read the properties of the physical bodies from the config file and write them into <Bodies>, a list of all physical objects of the simulation.
    Lightcurve = Bodies.calc_physics(P)  # Calculate body positions and the resulting lightcurve.
    animation = CurveSimAnimation(P, Bodies, Lightcurve)

# -*- coding: utf-8 -*-
# Curvesim - A Star System Lightcurve Simulator
# Curvesim calculates the movements and eclipses of celestial bodies and produces a video of this.
# Specify mass, radius and other properties of some stars and planets in a configuration file.
# Then run "curvesim.py <configfilename>" to produce the video.
# The video shows simultanously a view of the star system from the top and from the side and
# the lightcurve of the system's total luminosity over time.
# Usually you do not need to look at or even modify the python code. Instead control the program's
# outcome with the config file. The meaning of all program parameters is documented in the config file.
# Curvesim uses ffmpeg to convert the data into a video. Download ffmpeg from https://www.ffmpeg.org/download.html.
# Extract the zip file and add "<yourdriveandpath>\FFmpeg\bin" to Environment Variable PATH.
#
# Your questions and comments are welcome.
# Just open an issue on https://github.com/lichtgestalter/curvesim/issues to get my attention :)

from cs_animation import CurveSimAnimation
from cs_bodies import CurveSimBodies
from cs_parameters import CurveSimParameters

# If you run this script in a jupyter notebook, uncomment the next 2 lines in order to provide the name of your config file.
# import sys
# sys.argv[1] = "provide_your_config_file_name_here.ini"

def curvesim():
    parameters = CurveSimParameters()  # Read program parameters from config file.
    bodies = CurveSimBodies(parameters)  # Initialize the physical bodies, calculate their state vectors and generate their patches for the animation
    lightcurve = bodies.calc_physics(parameters)  # Calculate body positions and the resulting lightcurve.
    CurveSimAnimation(parameters, bodies, lightcurve)  # Create the video
    return parameters, bodies, lightcurve


if __name__ == '__main__':
    parameters_, bodies_, lightcurve_ = curvesim()
    # print(parameters_)
    print(bodies_)
    CurveSimBodies.prettyprint_lightcurve(lightcurve_)

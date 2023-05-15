# -*- coding: utf-8 -*-
# SSLS - Star System Lightcurve Simulator
# The SSLS calculates the movements and eclipses of celestial bodies and produces a video of this.
# Specify mass, radius and other properties of some stars and planets in a configuration file.
# Then run "cs_main.py <configfilename>" to produce the video.
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

# If you run this script in a jupyter notebook, uncomment this line in order to provide the name of your config file.
# sys.argv[1]="ssls.ini"

def main():
    config = configparser.ConfigParser(inline_comment_prefixes='#')
    configfilename = CurveSimParameters.find_and_check_config_file(default="ssls.ini")
    config.read(configfilename)
    p = CurveSimParameters(config, Standard_sections)  # Read program parameters from config file.
    bodies = CurveSimBodies(p, configfilename, Standard_sections)  # Initialize the physical bodies, calculate their state vectors and generate their patches for the animation
    lightcurve = bodies.calc_physics(p)  # Calculate body positions and the resulting lightcurve.
    CurveSimAnimation(p, bodies, lightcurve)  # Create the video


if __name__ == '__main__':
    main()

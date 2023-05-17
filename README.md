# curvesim
A Star System and Lightcurve Simulator<br>
<br>
Curvesim calculates the movements and eclipses of celestial bodies and produces a video of this and of the resulting lightcurve.<br>
<br>
Specify mass, radius and other properties of some stars and planets in a configuration file.<br>
Then run "curvesim.py <configfilename>" to produce the video.
The video shows simultanously a view of the star system from the top and from the side and
the lightcurve of the system's total luminosity over time.<br>
<br>
Usually you do not need to look at or even modify the python code. Instead control the program's
outcome with the config file. The meaning of all program parameters is documented in the config file.<br>
<br>
Curvesim uses ffmpeg to convert the data into a video. Download ffmpeg from https://www.ffmpeg.org/download.html.<br>
Extract the zip file and add "<yourdriveandpath>\FFmpeg\bin" to Environment Variable PATH.<br>
<br>
Your questions and comments are welcome.<br>
Just open an issue on https://github.com/lichtgestalter/curvesim/issues to get my attention :)<br>

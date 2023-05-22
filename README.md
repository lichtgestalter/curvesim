# curvesim
A Star System and Lightcurve Simulator<br>
<br>
Curvesim produces a video of the movements and eclipses of celestial bodies and of the resulting lightcurve.<br>
<br>
Curvesim is fast and the videos use little disk space. A video takes about the same time to produce as it is long and uses less than 0.5 MB disc space per minute.<br>
<br>
Specify mass, radius, orbital elements and other properties of some stars and planets in a configuration file.<br>
Then run "curvesim.py <configfilename>" to produce the video.
The video shows simultanously a view of the star system from the top and from the side and
the lightcurve of the system's total luminosity over time.<br>
<br>
Usually you do not need to look at or even modify the python code. Instead control the program's
outcome with the config file. The meaning of all program parameters is documented in the config file.<br>
<br>
Curvesim uses ffmpeg to convert the data into a video. <br> 
Download ffmpeg from https://www.ffmpeg.org/download.html. <br>
Extract the zip file and (on Windows) add "<yourdriveandpath>\FFmpeg\bin" to the environment variable PATH.<br>
<br>
For questions and comments just open an issue on https://github.com/lichtgestalter/curvesim/issues to get my attention :)<br>

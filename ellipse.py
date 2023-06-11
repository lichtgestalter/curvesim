import math
import numpy as np

def keplerian_elements_to_start_position(a, e, i, Ω, ω, nu):
    ea = 2 * math.atan(math.sqrt((1 - e) / (1 + e)) * math.tan(nu / 2))
    r = a * (1 - e * math.cos(ea))
    x = r * (math.cos(Ω) * math.cos(ω + nu) - math.sin(Ω) * math.sin(ω + nu) * math.cos(i))
    y = r * (math.sin(Ω) * math.cos(ω + nu) + math.cos(Ω) * math.sin(ω + nu) * math.cos(i))
    z = r * (math.sin(i) * math.sin(ω + nu))
    return np.array([x, y, z])




ea = 2 * atand(sqrt((1 - ee) / (1 + ee)) * tan(nu / 2))
r = a * (1 - ee * cos(ea))
xx = r * (cos(OM) * cos(om + nu) - sin(OM) * sin(om + nu) * cos(ii))
yy = r * (sin(OM) * cos(om + nu) + cos(OM) * sin(om + nu) * cos(ii))
zz = r * (sin(ii) * sin(om + nu))



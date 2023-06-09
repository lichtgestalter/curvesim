import numpy as np
from scipy.optimize import minimize

def fit_ellipse(points):
    def ellipse_distance(params):
        A, B, C, D, E, F = params
        distance = 0
        for i, (x, y, z) in enumerate(points):
            distance += ((A*x + B*y + C*z + D)**2) / (A**2 + B**2 + C**2) + ((E*x + F*y + 1 - z)**2)
        return distance

    def calculate_ellipse_points(params):
        A, B, C, D, E, F = params
        t = np.linspace(0, 2 * np.pi, 5)
        x_ellipse = (B*F - C*E) * np.sin(t) / (A*C - B**2) - (A*F - B*E) * np.cos(t) / (A*C - B**2)
        y_ellipse = (C*D - A*F) * np.sin(t) / (A*C - B**2) - (A*D - B*F) * np.cos(t) / (A*C - B**2)
        z_ellipse = (A - x_ellipse * E - y_ellipse * F) / C
        return list(zip(x_ellipse, y_ellipse, z_ellipse))

    initial_guess = (1, 1, 1, 1, 1, 1)  # Initial guess for ellipse parameters
    result = minimize(ellipse_distance, initial_guess, method='Nelder-Mead')

    parameters = result.x
    ellipse_points = calculate_ellipse_points(parameters)

    return parameters, ellipse_points



L0 = (51.16, 29.54, -34.11)
L72 = (-110.82, -63.98, 250.36)
L144 = (-252.55, -145.81, 236.72)
L216 = (-303.82, -175.41, 134.19)
L288 = (-243.18, -140.40, -14.36)

parameters, ellipse_points = fit_ellipse([L0, L72, L144, L216, L288])
print(parameters)
print(ellipse_points)

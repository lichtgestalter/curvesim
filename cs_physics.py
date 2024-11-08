import math


class CurveSimPhysics:

    @staticmethod
    def kepler_equation(ea, e, ma):
        """ea: eccentric anomaly [rad], e: eccentricity, ma: mean anomaly [rad]"""
        if not -2 * math.pi < ea < 2 * math.pi:
            raise ValueError("eccentric anomaly ea must be in radians but is outside of the range ]-2π;2π[")
        if not -2 * math.pi < ma < 2 * math.pi:
            raise ValueError("mean anomaly ma must be in radians but is outside of the range ]-2π;2π[")
        if not 0 <= e < 1:
            raise ValueError("eccentricity e is outside of the range [0;1[")
        return ea - e * math.sin(ea) - ma

    @staticmethod
    def kepler_equation_derivative(ea, e):
        """ea: eccentric anomaly [rad], e: eccentricity"""
        return 1.0 - e * math.cos(ea)

    @staticmethod
    def kepler_equation_root(e, ma, ea_guess=0.0, tolerance=1e-10, max_steps=50):
        """Calculate the root of the Kepler Equation with the Newton–Raphson method.
            e: eccentricity, ma: mean anomaly [rad], ea_guess: eccentric anomaly [rad]. ea_guess=ma is a good start."""
        for n in range(max_steps):
            delta = CurveSimPhysics.kepler_equation(ea_guess, e, ma) / CurveSimPhysics.kepler_equation_derivative(ea_guess, e)
            if abs(delta) < tolerance:
                return ea_guess - delta
            ea_guess -= delta
        raise RuntimeError('Newton\'s root solver did not converge.')

    @staticmethod
    def kepler_equation_root_chatgpt(e, ma, ea_guess=0.0, tolerance=1e-10, max_steps=50):
        """Returns the root of the Kepler Equation with the Newton–Raphson method.
            e: eccentricity, ma: mean anomaly [rad], ea_guess: eccentric anomaly [rad]. ea_guess=ma is a good start."""
        ea = ea_guess
        for _ in range(max_steps):
            f_ea = ea - e * math.sin(ea) - ma
            f_prime_ea = 1 - e * math.cos(ea)
            delta = -f_ea / f_prime_ea
            ea += delta
            if abs(delta) < tolerance:
                return ea
        raise RuntimeError("Kepler's equation solver did not converge.")

    @staticmethod
    def kepler_equation_root_debug(e, ma, ea_guess=0.0, tolerance=1e-10, max_steps=50):
        """
        Alternative Method for calculating the root of the Kepler Equation from source
        [f]: https://www.researchgate.net/publication/232203657_Orbital_Ephemerides_of_the_Sun_Moon_and_Planets, Section 8.10
        e: eccentricity, ma: mean anomaly [rad], ea_guess: eccentric anomaly [rad]. ea_guess=ma is a good start.
        """
        e_deg = math.degrees(e)
        ma_deg = math.degrees(ma)
        ea_deg = ma_deg + e_deg * math.sin(ma)

        for n in range(max_steps):
            delta_ma = ma_deg - (ea_deg - e * math.sin(ea_deg))
            delta_ea = delta_ma / (1 - e * math.cos(ea_deg))
            ea_deg += delta_ea
            if abs(delta_ea) < tolerance:
                return math.radians(ea_deg)
        raise RuntimeError('Solution for Kepler\'s Equation did not converge.')


    @staticmethod
    def gravitational_parameter(bodies, g):
        """Calculate the gravitational parameter of masses orbiting a common barycenter
        https://en.wikipedia.org/wiki/Standard_gravitational_parameter"""
        mass = 0.0
        for body in bodies:
            mass += body.mass
        return g * mass

    @staticmethod
    def distance_2d_ecl(body1, body2, i):
        """Return distance of the centers of 2 physical bodies as seen by a viewer (projection y->0)."""
        dx = body1.positions[i][0] - body2.positions[i][0]
        dz = body1.positions[i][2] - body2.positions[i][2]
        return math.sqrt((dx ** 2 + dz ** 2))

    @staticmethod
    def limbdarkening(relative_radius, beta):
        """https://en.wikipedia.org/wiki/Limb_darkening
        https://de.wikipedia.org/wiki/Photosph%C3%A4re#Mitte-Rand-Verdunkelung
        Approximates the flux of a star at a point on the star seen from a very large distance.
        The point's apparent distance from the star's center is relative_radius * radius.
        Beta depends on the wavelength. Beta=2.3 is a good compromise for the spectrum of visible light."""
        if relative_radius >= 1:
            return 1 / (1 + beta)
        return (1 + beta * math.sqrt(1 - relative_radius ** 2)) / (1 + beta)

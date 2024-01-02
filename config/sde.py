"""
Calculate input parameters for SDE
"""

import math
import plasmapy
import astropy.units as u
from astropy.constants.si import m_p, m_e, c, e
from plasmapy.particles import Particle
from plasmapy.formulary.frequencies import gyrofrequency
from plasmapy.formulary.lengths import gyroradius
from plasmapy.formulary.speeds import Alfven_speed, thermal_speed
from plasmapy.formulary.relativity import Lorentz_factor


def calc_sde_coefficients(eth,
                          bgauss,
                          lcorr,
                          knorm,
                          species,
                          va=1.0,
                          L0=1.0,
                          sigma2=1.0,
                          v0=1.0,
                          gamma=5 / 3):
    """Calculate diffusion coefficients

    Args:
        eth: thermal energy in keV
        bgauss: magnetic field in Gauss
        lcorr: correlation length in km
        knorm: kappa normalization in m^2/s
        species: particle species
        va: Alfven speed in m/s
        L0: length scale in m
        sigma2: turbulence variance
        v0: velocity normalization in km/s
        gamma: turbulence spectral slope
    """
    slen = 37
    print("Particle species:".ljust(slen), Particle(species))
    print("Particle thermal energy:".ljust(slen), eth)
    vth = thermal_speed(eth, species)
    print("Thermal speed:".ljust(slen), vth)
    print("Thermal speed/va:".ljust(slen), (vth / va))
    print("Thermal speed/v0:".ljust(slen), (vth / v0))
    eva = 0.5 * species.mass * va**2
    print("Particle energy with Alfven speed:".ljust(slen), eva.to("keV"))
    omega = gyrofrequency(bgauss,
                          species).to(u.s**-1,
                                      equivalencies=u.dimensionless_angles())
    print("Particle Gyrofrequency:".ljust(slen), omega)
    rg = gyroradius(bgauss, species, T=eth)
    print("Particle Gyroradius:".ljust(slen), rg)
    A0 = math.sin(math.pi / gamma) / (math.pi / gamma)
    kpara = vth**3 / (4 * math.pi * A0 * lcorr.to("m") * omega**2 * sigma2)
    kpara *= 1.0 + (omega * lcorr.to("m") / vth)**gamma * 8 / ((2 - gamma) *
                                                               (4 - gamma))
    print("kappa parallel:".ljust(slen), kpara)
    print("Normed kappa parallel:".ljust(slen), (kpara / knorm))
    duu0 = ((math.pi / 4) * A0 * sigma2 * omega**(2 - gamma) *
            (lcorr.to("m") / vth)**(1 - gamma))
    print("Duu0:".ljust(slen), duu0)
    print("Duu0*t0:".ljust(slen), duu0 * (L0 / v0))
    gamma_l = Lorentz_factor(vth)
    p = gamma_l * species.mass * vth
    v1 = e * v0 / (p * c)
    v2 = species.mass * e * v0 / p**2
    v1 *= bgauss.to("T") * L0
    v2 *= bgauss.to("T") * L0
    v1 = v1.to("", equivalencies=[u.C * u.T, u.kg * u.s**-1])
    v2 = v2.to("", equivalencies=[u.C * u.T, u.kg * u.s**-1])
    print("Parameter 1 for particle drift:".ljust(slen), v1)
    print("Parameter 2 for particle drift:".ljust(slen), v2)
    tau0 = 3 * kpara / vth**2  # scattering time scale
    tau0_scattering = tau0 / (L0 / va)
    print("Mean-free path:".ljust(slen), (tau0 * vth).to("km"))
    print("rg*vth/L0:".ljust(slen), (rg * vth / L0))
    print("Normed scattering time for particles:".ljust(slen), tau0_scattering)
    print("Normed mean free path for particles:".ljust(slen),
          (tau0 * vth / L0))
    print("omega*tau:".ljust(slen), (tau0 * omega))
    print("-" * slen)


def reconnection_test():
    """parameters for reconnection test
    """
    L0 = 5.0E6 * u.m
    b0 = 50.0 * u.G
    n0 = 1.0e10 * u.cm**-3
    gamma = 5.0 / 3.0
    lscale = 1.0
    bscale = 1.0
    slab_component = 1.0
    L0 *= lscale
    b0 *= bscale
    n0 *= bscale**2  # to keep the Alfven speed unchanged
    rho = n0 * (m_p + m_e)
    va = Alfven_speed(B=b0, density=rho)
    v0 = va
    knorm = L0 * v0
    Lc = 333 * u.km
    sigma2 = slab_component
    slen = 37  # string length
    print("Magnetic field normalization:".ljust(slen), b0)
    print("Number density normalization:".ljust(slen), n0)
    print("Length normalization:".ljust(slen), L0)
    print("Velocity normalization:".ljust(slen), v0)
    print("Time normalization:".ljust(slen), (L0 / v0))
    print("Alfven speed:".ljust(slen), va)
    print("Correlation length:".ljust(slen), Lc)
    print("Slab component:".ljust(slen), slab_component)
    print("Turbulence variance:".ljust(slen), sigma2)
    print("kappa normalization:".ljust(slen), knorm)
    print("Turbulence spectrum slope:".ljust(slen), gamma)
    print("-" * slen)
    calc_sde_coefficients(1 * u.keV, b0, Lc, knorm, Particle('e-'), va,
                          L0, sigma2, v0, gamma)


if __name__ == "__main__":
    reconnection_test()

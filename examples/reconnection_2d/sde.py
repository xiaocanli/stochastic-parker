"""
Analysis functions for SDE calculation.
We use SI units here
"""

import math

UNIT_CHARGRE = 1.6021765E-19
LIGHT_SPEED = 299792458.0


def calc_ptl_speed(eth, species, verbose=False):
    """ calculate thermal speed

    Args:
        eth: thermal energy in keV
        species: particle species
    """
    if species == 'e':
        pmass = 9.10938356E-31  # in kilogram
        ptl_name = "electron"
    else:
        pmass = 1.6726219E-27  # in kilogram
        ptl_name = "proton"
    rest_energy = pmass * LIGHT_SPEED**2
    gamma = eth * 1E3 * UNIT_CHARGRE / rest_energy + 1
    vth = math.sqrt(1.0 - 1.0 / gamma**2)  # in m/s
    if verbose:
        print("Thermal speed for %0.1f keV %s: %fc" % (eth, ptl_name, vth))
    return vth


def calc_va_si(b0, n0, verbose=False):
    """Calculate the Alfven speed

    Args:
        b0: magnetic field strength in Gauss
        n0: particle number density in cm^-3
    """
    pmass = 1.6726219E-27  # in kilogram
    mu0 = 4 * math.pi * 1E-7
    va = b0 * 1E-4 / math.sqrt(mu0 * n0 * 1E6 * pmass)
    print("The Alfven speed is %f" % va)
    return va


def calc_va_cgs(n, B):
    """Calculate the Alfven speed

    Args:
        n: density (cm^-3)
        B: magnetic field (Gauss)
    """
    mi = 1.6726219E-24  # proton mass in gram
    va = B / math.sqrt(4 * math.pi * n * mi)  # CGS in cm/s
    return va


def calc_gyro_frequency(bgauss, species, verbose=False):
    """Calculate particle gyro-frequency

    Args:
        bgauss: magnetic field in Gauss
        species: particle species
    """
    if species == 'e':
        pmass = 9.10938356E-31  # in kilogram
        ptl_name = "electron"
    else:
        pmass = 1.6726219E-27  # in kilogram
        ptl_name = "proton"
    omega = UNIT_CHARGRE * bgauss * 1E-4 / pmass
    if verbose:
        print("The gyro-frequency of %s in %0.1f Gauss magnetic field: %e Hz" %
              (ptl_name, bgauss, omega))
    return omega


def calc_kappa_parallel(eth,
                        bgauss,
                        lcorr,
                        knorm,
                        species,
                        va=1.0,
                        L0=1.0,
                        sigma2=1.0):
    """Calculate parallel diffusion coefficient

    Args:
        eth: thermal energy in keV
        bgauss: magnetic field in Gauss
        lcorr: correlation length in km
        knorm: kappa normalization in m^2/s
        species: particle species
        va: Alfven speed in m/s
        L0: length scale in m
        sigma2: turbulence variance
    """
    csc_3pi_5 = 1.05146
    lcorr *= 1000  # in meter
    if species == 'e':
        pmass = 9.10938356E-31  # in kilogram
        ptl_name = "electron"
    else:
        pmass = 1.6726219E-27  # in kilogram
        ptl_name = "proton"
    vth = calc_ptl_speed(eth, species) * LIGHT_SPEED  # in m/s
    print("Magnetic field: %0.1f G" % bgauss)
    print(ptl_name + " energy: %0.1f keV" % eth)
    print("Thermal speed: %e km/s" % (vth / 1E3))
    print("Thermal speed/va: %f" % (vth / va))
    eva = 0.5 * pmass * va**2 / UNIT_CHARGRE / 1000  # in kev
    print("Particle energy with speed va: %f keV" % (eva))
    omega = calc_gyro_frequency(bgauss, species)
    rg = vth / omega  # gyro-radius in m
    kpara = 3 * vth**3 * csc_3pi_5 / (20 * lcorr * omega**2 * sigma2)
    kpara *= 1.0 + (72.0 / 7) * (omega * lcorr / vth)**(5 / 3)
    # kpara = 2 * vth**3 / (3 * omega**2 * math.pi)
    kperp = (5 * vth * lcorr / 12) * math.sin(3 * math.pi / 5) * sigma2
    gamma = 1.0 / math.sqrt(1 - (vth / LIGHT_SPEED)**2)
    p = gamma * pmass * vth
    v1 = UNIT_CHARGRE * va / (p * LIGHT_SPEED)
    v2 = pmass * UNIT_CHARGRE * va / p**2
    v1 *= bgauss * 1E-4 * L0
    v2 *= bgauss * 1E-4 * L0
    tau0 = 3 * kpara / vth**2  # scattering time scale
    tau0_scattering = tau0 / (L0 / va)
    print("Gyro-frequency : %e Hz" % omega)
    print("Gyro-radius: %e m" % rg)
    print("Mean-free path: %e m" % (tau0 * vth))
    print("kappa parallel: %e m^2/s" % kpara)
    print("Normed kappa parallel: %e" % (kpara / knorm))
    print("kperp / kpara: %e" % (kperp / kpara))
    print("Parameters for particle drift: %e, %e" % (v1, v2))
    print("rg*vth/L0: %e" % (rg * vth / L0))
    print("Scattering time for initial particles: %e" % tau0_scattering)
    print("Mean free path initial particles: %e\n" % (tau0 * vth))


def reconnection_test():
    """parameters for reconnection test
    """
    L0 = 5.0E6  # in m
    b0 = 50.0  # in Gauss
    n0 = 1.0e10  # number density in cm^-3
    lscale = 1.0  # scale the MHD length
    bscale = 1.0  # scale the magnetic field
    slab_component = 1.0
    L0 *= lscale
    b0 *= bscale
    n0 *= bscale**2  # to keep the Alfven speed unchanged
    va = calc_va_cgs(n0, b0) / 1E2  # in m/s
    knorm = L0 * va
    Lc = 333  # in km
    sigma2 = slab_component
    print("Alfven speed: %0.1f km/s" % (va / 1E3))
    print("Correlation length: %0.1f km" % Lc)
    print("Turbulence variance: %e\n" % sigma2)
    calc_kappa_parallel(1, b0, Lc, knorm, 'e', va, L0, sigma2)


if __name__ == "__main__":
    reconnection_test()

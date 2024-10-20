The Equation Solved
-------------------

Parker’s transport equation

.. math::
   :name: equ_parker

   \frac{\partial f}{\partial t} + (\boldsymbol{V}+\boldsymbol{V}_d)\cdot\nabla f
     - \frac{1}{3}\nabla\cdot\boldsymbol{V}\frac{\partial f}{\partial\ln p}
     = \nabla\cdot(\boldsymbol{\kappa}\nabla f) + Q,

where :math:`f(x_i, p, t)` is the particle distribution function as a function of the particle position :math:`x_i`,  momentum :math:`p` (isotropic momentum assumed), and time :math:`t`; :math:`\boldsymbol{\kappa}` is the spatial diffusion coefficient tensor, :math:`\boldsymbol{V}` is the bulk plasma velocity, :math:`\boldsymbol{V}_d=\frac{pcw}{3q}\nabla\times\left(\frac{\boldsymbol{B}}{B^2}\right)` is the particle drift, where :math:`p=\gamma m v` is particle momentum, :math:`c` is the speed of light, :math:`w=v/c` is the normalized particle speed, and :math:`q` is particle charge.

.. note::

    The particle drift velocity is only needed in the 3D model since it is out-of-plane in a 2D model.

and :math:`Q` is the source. In general, :math:`\boldsymbol{V}` can be obtained directly from the MHD simulations, :math:`\boldsymbol{V}_d` and :math:`\boldsymbol{\kappa}` both depend on the vector magnetic field, and :math:`Q` could depend on many plasma properties (e.g., number density, current density, and compression). :math:`\boldsymbol{\kappa}` also depends on the turbulence properties (e.g., amplitude, spectral slope, and anisotropy). The diffusion coefficient tensor is given by

.. math::

   \kappa_{ij} = \kappa_\perp\delta_{ij} -
     \frac{(\kappa_\perp-\kappa_\parallel)B_iB_j}{B^2},

where :math:`\kappa_\parallel` and :math:`\kappa_\perp` are the parallel
and perpendicular diffusion coefficients. Here :math:`\kappa_\parallel`
can be calculated from the quasi-linear
theory [Jokipii71]_. Assuming that
magnetic turbulence is well developed and has an isotropic power
spectrum :math:`P\sim k^{-5/3}`, the resulting
:math:`\kappa_\parallel\sim p^{4/3}` when the particle gyroradius is
much smaller than the correlation length of turbulence. In particular,
we use the following expression for
:math:`\kappa_\parallel` [Giacalone99]_,

.. math::
   :name: equ_kpara_qlt

   \begin{aligned}
     \kappa_\parallel(v) & = \frac{3v^3}{20L_c\Omega_0^2\sigma^2}
     \csc\left(\frac{3\pi}{5}\right)\left[1+\frac{72}{7}
     \left(\frac{\Omega_0L_c}{v}\right)^{5/3}\right]\\
     & \approx1.622\frac{v^{4/3}L_c^{2/3}}{\Omega_0^{1/3}\sigma^2}
   \end{aligned}

where :math:`v` is the particle speed, :math:`L_c` is the turbulence
correlation length, :math:`\Omega_0` is the particle gyrofrequency, and
:math:`\sigma^2=\left<\delta B^2\right>/B_0^2` is the normalized wave
variance of turbulence. Reference [Giacalone99]_ gave
a derivation of equ_kpara_qlt_ . Below is
summary of it with some missing details. The velocity-dependent parallel
diffusion coefficient is

.. math:: \kappa_\parallel(v) = \frac{v^2}{4}\int_0^1 \frac{(1-\mu^2)^2d\mu}{D_{\mu\mu}}

where :math:`\mu` is cosine of the pitch angle and :math:`D_{\mu\mu}`
the pitch-angle diffusion coefficient. :math:`D_{\mu\mu}` is related to
magnetic field fluctuations.

.. math:: D_{\mu\mu} = \frac{\pi}{4}\Omega_0(1-\mu^2)\frac{k_\text{res}P(k_\text{res})}{B_0^2}

where :math:`\Omega_0` is particle gyrofrequency,
:math:`k_\text{res}=|\Omega_0/v\mu|` is the resonant wavenumber. The
above equation is strictly applicable only for the case of 1D turbulence
in which the wavevectors are aligned with the mean field. For
anisotropic turbulence (e.g., 2D + slab), only the slab turbulence
(about 20% of all turbulent
fluctuations [Bieber96]_) affect particle
parallel transport [Florinski03]_. The
turbulence power is usually expressed as

.. math:: P(k) = \frac{\left<\delta B^2\right>}{1+(kL_c)^\gamma}\left[\int_{k_\text{min}}^{k_\text{max}}\frac{dk}{1+(kL_c)^\gamma}\right]^{-1}

where :math:`k_\text{min}` and :math:`k_\text{max}` are the smallest and
largest wavenumbers in the system, :math:`L_c` is the turbulence
correlation length, and :math:`\gamma` is the turbulence spectrum index
(e.g., 5/3). For :math:`k_\text{min}\ll 1/L_c \ll k_\text{max}`, the
integral in the above equation can be taken from 0 to :math:`\infty`.
From the table of integral,
:math:`\int_0^\infty x^{\mu-1} dx / (1+x^\nu) = \pi\csc(\mu\pi/\nu)/\nu`.
Then,

.. math:: P(k) = \frac{\left<\delta B^2\right>L_c}{1+(kL_c)^\gamma}\left[(\pi/\gamma)\csc(\pi/\gamma)\right]^{-1}

The parallel diffusion coefficient is

.. math::

   \begin{aligned}
     \kappa_\parallel(v) & = \frac{v^2\csc(\pi/\gamma)}{\Omega_0^2\sigma^2\gamma L_c}\int_0^1(1-\mu^2)|v\mu|\left(1+\left|\frac{\Omega_0}{v\mu}\right|L_c\right)^\gamma d\mu \\
     & = \frac{v^3\csc(\pi/\gamma)}{4\Omega_0^2\sigma^2\gamma L_c}\left[1+ \left(\frac{\Omega_0L_c}{v}\right)^\gamma\frac{8}{(2-\gamma)(4-\gamma)}\right]
   \end{aligned}

where we assume :math:`v>0` and :math:`\gamma<2`. When
:math:`\gamma=5/3` (Kolmogorov), we can get
equ_kpara_qlt_. (What will happen when
:math:`\gamma>2`?).

Test-particle simulations have suggested that
:math:`\kappa_\perp/\kappa_\parallel` is about 0.02-0.04 and is nearly
independent of particle
energy [Giacalone99]_. There is also
observational evidence suggesting that
:math:`\kappa_\perp/\kappa_\parallel` can be quite
large [Dwyer97]_ [Zhang03]_.

The Parker transport equation can be solved by integrating the
stochastic differential equation corresponding to the Fokker–Planck form
of the transport
equation [Zhang99]_ [Florinski09]_ [Pei10]_ [Kong17]_.
Neglecting the source term :math:`Q` in
equ_parker_ and assuming :math:`F=fp^2`,

.. math::

   \begin{aligned}
     \frac{\partial F}{\partial t}
     & = -\nabla\cdot\left[(\nabla\cdot\boldsymbol{\kappa}+\boldsymbol{V})F\right] +
     \frac{\partial}{\partial p} \left[\frac{p}{3}\nabla\cdot\boldsymbol{V} F\right] +
     \nabla\cdot(\nabla\cdot(\boldsymbol{\kappa}F)),
   \end{aligned}

which is equivalent to a system of stochastic differential equations
(SDEs) of the Ito type.


.. [Bieber96] Bieber, J.W., Wanner, W. and Matthaeus, W.H., 1996. Dominant two‐dimensional solar wind turbulence with implications for cosmic ray transport. Journal of Geophysical Research: Space Physics, 101(A2), pp.2511-2522.
.. [Dwyer97] Dwyer, J.R., Mason, G.M., Mazur, J.E., Jokipii, J.R., Von Rosenvinge, T.T. and Lepping, R.P., 1997. Perpendicular transport of low-energy corotating interaction region-associated nuclei. The Astrophysical Journal, 490(1), p.L115.
.. [Florinski03] Florinski, V., Zank, G.P. and Pogorelov, N.V., 2003. Galactic cosmic ray transport in the global heliosphere. Journal of Geophysical Research: Space Physics, 108(A6).
.. [Florinski09] Florinski, V. and Pogorelov, N.V., 2009. Four-dimensional transport of galactic cosmic rays in the outer heliosphere and heliosheath. The Astrophysical Journal, 701(1), p.642.
.. [Giacalone99] Giacalone, J. and Jokipii, J.R., 1999. The transport of cosmic rays across a turbulent magnetic field. The Astrophysical Journal, 520(1), p.204.
.. [Jokipii71] Jokipii, J.R., 1971. Propagation of cosmic rays in the solar wind. Reviews of Geophysics, 9(1), pp.27-87.
.. [Kong17] Kong, X., Guo, F., Giacalone, J., Li, H. and Chen, Y., 2017. The acceleration of high-energy protons at coronal shocks: the effect of large-scale streamer-like magnetic field structures. The Astrophysical Journal, 851(1), p.38.
.. [Pei10] Pei, C., Bieber, J. W., Burger, R. A., & Clem, J. 2010, Journal of Geophysical Research (Space Physics), 115, A12107
.. [Zhang99] Zhang, M., 1999. A Markov stochastic process theory of cosmic-ray modulation. The Astrophysical Journal, 513(1), p.409.
.. [Zhang03] Zhang, M., Jokipii, J.R. and McKibben, R.B., 2003. Perpendicular transport of solar energetic particles in heliospheric magnetic fields. The Astrophysical Journal, 595(1), p.493.
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

.. math:: P(k) = \frac{\left<\delta B^2\right>_\text{slab}}{1+(kl_\text{slab})^\gamma}\left[\int_{k_\text{min}}^{k_\text{max}}\frac{dk}{1+(kl_\text{slab})^\gamma}\right]^{-1}

where :math:`k_\text{min}` and :math:`k_\text{max}` are the smallest and
largest wavenumbers in the system, :math:`l_\text{slab}` is the turbulence
correlation length for the slab component, and :math:`\gamma` is the turbulence spectrum index
(e.g., 5/3). For :math:`k_\text{min}\ll 1/l_\text{slab} \ll k_\text{max}`, the
integral in the above equation can be taken from 0 to :math:`\infty`.
From the table of integral,
:math:`\int_0^\infty x^{\mu-1} dx / (1+x^\nu) = \pi\csc(\mu\pi/\nu)/\nu`.
Then,

.. math:: P(k) = \frac{\left<\delta B^2\right>l_\text{slab}}{1+(kl_\text{slab})^\gamma}\left[(\pi/\gamma)\csc(\pi/\gamma)\right]^{-1}

The parallel diffusion coefficient is

.. math::

   \begin{aligned}
     \kappa_\parallel(v) & = \frac{v^2\csc(\pi/\gamma)}{\Omega_0^2\sigma^2\gamma l_\text{slab}}\int_0^1(1-\mu^2)|v\mu|\left(1+\left|\frac{\Omega_0}{v\mu}\right|l_\text{slab}\right)^\gamma d\mu \\
     & = \frac{v^3\csc(\pi/\gamma)}{4\Omega_0^2\sigma^2\gamma l_\text{slab}}\left[1+ \left(\frac{\Omega_0l_\text{slab}}{v}\right)^\gamma\frac{8}{(2-\gamma)(4-\gamma)}\right]
   \end{aligned}

where we assume :math:`v>0` and :math:`\gamma<2`. When
:math:`\gamma=5/3` (Kolmogorov), we can get
equ_kpara_qlt_. (What will happen when
:math:`\gamma>2`?).

The resulting parallel mean free path is

.. math::

  \begin{align}
    \lambda_\parallel & =3\kappa_\parallel/v\\
    & =\frac{3v^2\csc(\pi/\gamma)}{4\Omega_0^2\sigma_\text{slab}^2\gamma l_\text{slab}}\left[1+ \left(\frac{\Omega_0l_\text{slab}}{v}\right)^\gamma\frac{8}{(2-\gamma)(4-\gamma)}\right]\\
    & = \frac{3\csc(\pi/\gamma)}{\gamma\sigma_\text{slab}^2}l_\text{slab}R^2\left[\frac{1}{4}+\frac{2R^{-\gamma}}{(2-\gamma)(4-\gamma)}\right]
  \end{align}

where :math:`R=R_L/l_\text{slab}` with :math:`R_L=v/\Omega_0` the Larmor radius for a :math:`90^\circ` pitch angle. Test-particle simulations have suggested that :math:`\kappa_\perp/\kappa_\parallel` is about 0.02-0.04 at 1AU and is nearly independent of particle energy [Giacalone99]_. There is also observational evidence suggesting that :math:`\kappa_\perp/\kappa_\parallel` can be quite large [Dwyer97]_ [Zhang03]_. Theoretically, :math:`\kappa_\perp` can be calculated based nonlinear guiding center theory (NLGC) [Matthaeus03]_ [Wijsen23]_.

.. math::

  \kappa_\perp=\mu^2v\lambda_\parallel^{1/3}\left(\alpha^2\sqrt{3\pi}\frac{2\nu-1}{\nu}\frac{\Gamma(\nu)}{\Gamma(\nu-1/2)}\sigma_\text{2D}^2l_\text{2D}\right)^{2/3}
  
  
where :math:`\nu=\gamma/2`, :math:`\Gamma` is the gamma function, and :math:`\alpha` is the free parameter, which is :math:`\sqrt{1/3}` according to [Matthaeus03]_. Then, the perpendicular mean free path is

.. math::

    \lambda_\perp=\frac{3}{2v}\int_{-1}^1\kappa_\perp d\mu
    = \lambda_\parallel^{1/3}\left(\alpha^2\sqrt{3\pi}\frac{2\nu-1}{\nu}\frac{\Gamma(\nu)}{\Gamma(\nu-1/2)}\sigma_\text{2D}^2l_\text{2D}\right)^{2/3}

For SEP transport, the unknowns left are the radial dependence of the turbulence properties (:math:`\sigma^2`, :math:`l_\text{slab}`, :math:`l_\text{2D}`). For :math:`\sigma^2`, we can prescribed it. For example, [Wijsen23]_ gave the following form.

.. math::

    \sigma^2=
        \begin{cases}
        \Lambda_0(r/r_0)^{\alpha_1} & r\leq r_1=0.5 \text{ au}\\
        \Lambda_1(r/r_0)^{\alpha_2} & r_1<r\leq r_2=2.0 \text{ au}\\
        \Lambda_2 & r>r_2
        \end{cases}

where :math:`r_0=0.1` au, :math:`\Lambda_0=0.1`, :math:`\Lambda_1=\Lambda_0(r_1/r_0)^{\alpha_1-\alpha_2}\approx0.15`, :math:`\Lambda_2=\Lambda_0(r_2/r_0)^{\alpha_2}\approx0.32`, with :math:`\alpha_1=0.5` and :math:`\alpha_2=0.25`. The choice of the above equation is to capture the different radial dependence of the background magnetic field and the turbulent field :math:`\delta B`. The 2D and slab correlation length are prescribed as :math:`l_\text{2D}=(0.0074 \text{ au}) (r/1 \text{ au})^{1.1}` and :math:`l_\text{slab}=3.9\times l_\text{2D}` in [Wijsen23]_.

For the transport modeling, :math:`\kappa_\parallel` part can be the same as before. The problem is :math:`\kappa_\perp`, which depends on both the 2D and the slab components. :math:`\kappa_\perp` also has different velocity dependence. Since we want to hide the analytical details, it is better not to include all the calculation of the diffusion coefficients in the code.

.. math::

    \begin{align}
    \kappa_\parallel &\sim v^{3-\gamma}B_0^{\gamma-2}(\sigma_\text{slab}^2)^{-1}l_\text{slab}^{\gamma-1} \\
    \lambda_\parallel &\sim v^{2-\gamma}B_0^{\gamma-2}(\sigma_\text{slab}^2)^{-1}l_\text{slab}^{\gamma-1} \\
    \kappa_\perp &\sim v^{(5-\gamma)/3}B_0^{(\gamma-2)/3}(\sigma_\text{slab}^2)^{-1/3}(\sigma_\text{2D}^2)^{2/3}l_\text{slab}^{(\gamma-1)/3}l_\text{2D}^{2/3}
    \end{align}

If we assume :math:`\sigma_\text{slab}^2 \sim \sigma_\text{2D}^2\sim\sigma^2` and :math:`l_\text{slab}\sim l_\text{2D}\sim L_c`,

.. math::

    \begin{align}
    \kappa_\parallel &\sim v^{3-\gamma}B_0^{\gamma-2}(\sigma^2)^{-1}L_c^{\gamma-1} \\
    \kappa_\perp &\sim v^{(5-\gamma)/3}B_0^{(\gamma-2)/3}(\sigma^2)^{1/3}L_c^{(\gamma+1)/3}
    \end{align}

The gradient of the diffusion coefficient can be calculated as

.. math::

    \frac{d\kappa_\parallel}{dx} \sim \kappa_\parallel\left(\frac{\gamma-2}{B_0}\frac{\partial B_0}{\partial x} - \frac{1}{\sigma_\text{slab}^2}\frac{\partial(\sigma_\text{slab}^2)}{\partial x} + (\gamma-1)\frac{1}{l_\text{slab}}\frac{\partial l_\text{slab}}{\partial x}\right)


.. math::

    \frac{d\kappa_\perp}{dx}\sim \kappa_\perp\left(\frac{\gamma-2}{3B_0}\frac{\partial B_0}{\partial x} - \frac{1}{3\sigma_\text{slab}^2}\frac{\partial(\sigma_\text{slab}^2)}{\partial x} + \frac{2}{3\sigma_\text{2D}^2}\frac{\partial(\sigma_\text{2D}^2)}{\partial x} + \frac{\gamma-1}{3}\frac{1}{l_\text{slab}}\frac{\partial l_\text{slab}}{\partial x} + \frac{2}{3}\frac{1}{l_\text{2D}}\frac{\partial l_\text{2D}}{\partial x}\right)


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
.. [Matthaeus03] Matthaeus, William H., Gang Qin, John William Bieber, and Gary Paul Zank. "Nonlinear collisionless perpendicular diffusion of charged particles." The Astrophysical Journal 590, no. 1 II (2003).
.. [Pei10] Pei, C., Bieber, J. W., Burger, R. A., & Clem, J. 2010, Journal of Geophysical Research (Space Physics), 115, A12107
.. [Wijsen23] Wijsen, Nicolas, Gang Li, Zheyi Ding, David Lario, Stefaan Poedts, Rachael Jo Filwett, Robert Colby Allen, and M. A. Dayeh. "On the seed population of solar energetic particles in the inner heliosphere." Journal of Geophysical Research: Space Physics 128, no. 3 (2023): e2022JA031203.
.. [Zhang99] Zhang, M., 1999. A Markov stochastic process theory of cosmic-ray modulation. The Astrophysical Journal, 513(1), p.409.
.. [Zhang03] Zhang, M., Jokipii, J.R. and McKibben, R.B., 2003. Perpendicular transport of solar energetic particles in heliospheric magnetic fields. The Astrophysical Journal, 595(1), p.493.

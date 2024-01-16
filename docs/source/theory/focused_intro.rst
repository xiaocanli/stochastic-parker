The Equation Solved
-------------------

.. warning:: 
  This part is still in early development. Please use it with caution.

.. note:: 
  It is better to read the section on Parker's transport equation first since some of the details are skipped here.

The focused transport equation is [Zank14]_

.. math::

   \begin{aligned}
     & \frac{\partial f}{\partial t} + (U_i + c\mu b_i)\frac{\partial f}{\partial x_i} + \frac{dc}{dt}\frac{\partial f}{\partial c} + \frac{d\mu}{dt}\frac{\partial f}{\partial\mu} = \left<\left.\frac{\delta f}{\delta t}\right\vert_s\right> \\
     & \frac{dc}{dt} = \left[\frac{1-3\mu^2}{2}b_ib_j\frac{\partial U_i}{\partial x_j}-\frac{1-\mu^2}{2}\nabla\cdot\boldsymbol{U}-\frac{\mu b_i}{c}\left(\frac{\partial U_i}{\partial t} + U_j\frac{\partial U_i}{\partial x_j}\right)\right]c \\
     & \frac{d\mu}{dt} = \frac{1-\mu^2}{2}\left[c\nabla\cdot\boldsymbol{b}+\mu\nabla\cdot\boldsymbol{U} - 3\mu b_ib_j\frac{\partial U_i}{\partial x_j}-\frac{2b_i}{c}\left(\frac{\partial U_i}{\partial t} + U_j\frac{\partial U_i}{\partial x_j}\right)\right]
   \end{aligned}

where :math:`\boldsymbol{c}` is the particle velocity in the flow frame
(:math:`\boldsymbol{v}=\boldsymbol{c} + \boldsymbol{U}`),
:math:`\mu\equiv\cos\theta=\boldsymbol{c}\cdot\boldsymbol{b}/c` is the
particle pitch angle, :math:`\boldsymbol{b}\equiv\boldsymbol{B}/B` is
the unit vector along the magnetic field. The following focused
transport equation is often
used [Zhang09]_ [Zuo13]_ [Zhang17]_ [Kong22]_,

.. math::
   :name: equ_ft

   \begin{aligned}
     \frac{\partial f}{\partial t} = \nabla\cdot(\boldsymbol{\kappa}_\perp\nabla f) - (v\mu\boldsymbol{b} + \boldsymbol{V} + \boldsymbol{V}_d)\cdot\nabla f + \frac{\partial}{\partial\mu}\left(D_{\mu\mu}\frac{\partial f}{\partial\mu}\right) - \frac{d\mu}{dt}\frac{\partial f}{\partial\mu} - \frac{dp}{dt}\frac{\partial f}{\partial p}
   \end{aligned}

where the terms on the right-hand side are cross-field spatial diffusion
with a tensor
:math:`\boldsymbol{\kappa}_\perp=\kappa_\perp\left(\overline{\overline{\mathbf{I}}}-\boldsymbol{b}\boldsymbol{b}\right)` [Zhang09]_,
streaming along the ambient or average magnetic field direction
:math:`\boldsymbol{b}` with particle speed :math:`v` and pitch-angle
cosine :math:`\mu`, convection with the background plasma
:math:`\boldsymbol{V}`, partial gradient/curvature drift
:math:`\boldsymbol{V}_d`, pitch-angle diffusion with a coefficient
:math:`D_{\mu\mu}`, focusing :math:`d\mu/dt`, and adiabatic
heating/cooling :math:`dp/dt`. In the adiabatic approximation, the drift
velocity, focusing rate, and cooling rate may be calculated from the
ambient magnetic field :math:`\boldsymbol{B}=B\boldsymbol{b}` and plasma
velocity :math:`\boldsymbol{V}` through

.. math::

   \begin{aligned}
     & \boldsymbol{V}_d=\frac{cpv}{qB}\left\{\frac{1-\mu^2}{2}\frac{\boldsymbol{B}\times\nabla B}{B^2}+\mu^2\frac{\boldsymbol{B}\times[(\boldsymbol{B}\cdot\nabla)\boldsymbol{B}]}{B^3}+\frac{1-\mu^2}{2}\frac{\boldsymbol{B}(\boldsymbol{B}\cdot\nabla\times\boldsymbol{B})}{B^3}\right\}\\
     & \frac{d\mu}{dt} = \frac{1-\mu^2}{2}\left[-v\boldsymbol{b}\cdot\nabla\ln B+\mu\nabla\cdot\boldsymbol{V} - 3\mu b_ib_j\frac{\partial V_i}{\partial x_j}-\frac{2b_i}{v}\frac{dV_i}{dt}\right] \\
     & \frac{dp}{dt} = -p\left[\frac{1-\mu^2}{2}\left(\nabla\cdot\boldsymbol{V}-b_ib_j\frac{\partial V_i}{\partial x_j}\right)+\mu^2b_ib_j\frac{\partial V_i}{\partial x_j}+\frac{\mu b_i}{v}\frac{dV_i}{dt}\right]
   \end{aligned}

where :math:`\boldsymbol{V}_d` includes gradient, curvature, and
parallel drifts. We may ignore :math:`\partial V_i/\partial t` in
:math:`dV_i/dt` if the flow is not dramatically evolving.

In the quasi-linear theory, resonant interaction between the particle
and the turbulent magnetic field can be related by the pitch-angle
diffusion coefficient
:math:`D_{\mu\mu}` [Jokipii77]_ [Kong22]_

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

.. math:: P(k) = \frac{\left<\delta B^2\right>L_c}{1+(kL_c)^\gamma}A_0

where :math:`A_0=\left[(\pi/\gamma)\csc(\pi/\gamma)\right]^{-1}`. In the
nonrelativistic limit, we take the pitch-angle diffusion coefficient in
the form of [Kong22]_

.. math::

   \begin{aligned}
     D_{\mu\mu}=D_{\mu\mu0}\left(\frac{p}{p_0}\right)^{\gamma-1}(1-\mu^2)(|\mu|^{\gamma-1} + h_0)
   \end{aligned}

where
:math:`D_{\mu\mu0}=\frac{\pi}{4}A_0\sigma^2\Omega_0^{2-\gamma}L_c^{1-\gamma}v_0^{\gamma-1}`,
and :math:`p_0`\ (:math:`v_0`) is the initial particle momentum
(velocity) at the injection energy. The parameter :math:`h_0` is added
to describe the scattering through :math:`\mu=0`, and we set
:math:`h_0=0.2` [Zhang17]_ [Kong22]_.

.. math::

   \begin{aligned}
     \frac{\partial D_{\mu\mu}}{\partial\mu} = D_{\mu\mu0}\left(\frac{p}{p_0}\right)^{\gamma-1}\left[-2\mu(|\mu|^{\gamma-1} + h_0) + (1-\mu^2)\text{sign}(\mu)|\mu|^{\gamma-2}\right]
   \end{aligned}

.. [Kong22] Kong, X., Chen, B., Guo, F., Shen, C., Li, X., Ye, J., Zhao, L., Jiang, Z., Yu, S., Chen, Y. and Giacalone, J., 2022. Numerical modeling of energetic electron acceleration, transport, and emission in solar flares: connecting loop-top and footpoint hard X-ray sources. The Astrophysical Journal Letters, 941(2), p.L22.
.. [Zank14] Zank, G.P., 2014. Transport processes in space physics and astrophysics (Vol. 877, p. 185). Berlin: Springer.
.. [Zhang09] Zhang, M., Qin, G. and Rassoul, H., 2009. Propagation of solar energetic particles in three-dimensional interplanetary magnetic fields. The Astrophysical Journal, 692(1), p.109.
.. [Zhang17] Zhang, M. and Zhao, L., 2017. Precipitation and release of solar energetic particles from the solar coronal magnetic field. The Astrophysical Journal, 846(2), p.107.
.. [Zuo13] Zuo, P., Zhang, M. and Rassoul, H.K., 2013. The role of cross-shock potential on pickup ion shock acceleration in the framework of focused transport theory. The Astrophysical Journal, 776(2), p.93.
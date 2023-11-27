Momentum Diffusion
------------------

.. note:: 
  Momentum diffusion is only partially supported now.

We can include an momentum diffusion term to the right side of
equ_parker_.

.. math::
   :name: equ_parker_2nd

   \frac{\partial f}{\partial t} + (\boldsymbol{V}+\boldsymbol{V}_d)\cdot\nabla f
     - \frac{1}{3}\nabla\cdot\boldsymbol{V}\frac{\partial f}{\partial\ln p}
     = \nabla\cdot(\boldsymbol{\kappa}\nabla f) +
     \frac{1}{p^2}\frac{\partial}{\partial p}
     \left(p^2D_{pp}\frac{\partial f}{\partial p}\right) + Q,

Neglecting the source term :math:`Q` in
equ_parker_2nd_ and assuming :math:`F=fp^2`,

.. math::

   \begin{aligned}
     \frac{\partial F}{\partial t} =
     & -\nabla\cdot\left[(\nabla\cdot\boldsymbol{\kappa}+\boldsymbol{V}+\boldsymbol{V}_d)F\right] +
     \nabla\cdot(\nabla\cdot(\boldsymbol{\kappa}F)) + \nonumber \\
     & \frac{\partial}{\partial p} \left[\left(\frac{p}{3}\nabla\cdot\boldsymbol{V} -
     \frac{\partial D_{pp}}{\partial p} - \frac{2D_{pp}}{p}\right) F\right] +
     \frac{\partial(D_{pp}F)}{\partial p^2}.
   \end{aligned}

which is equivalent to a system of SDEs of the Ito type,

.. math::

   \begin{aligned}
     dX & = (\nabla\cdot\boldsymbol{\kappa} + \boldsymbol{V} + \boldsymbol{V}_d)ds +
     \sum_\sigma\boldsymbol{\alpha}_\sigma dW_\sigma(s) \\
     dp & = \left(-\frac{p}{3}\nabla\cdot\boldsymbol{V} +
     \frac{\partial D_{pp}}{\partial p} + \frac{2D_{pp}}{p}\right)ds +
     \sqrt{2D_{pp}}dW(s)
   \end{aligned}

where
:math:`\sum_\sigma\alpha_\sigma^\mu\alpha_\sigma^\nu = 2\kappa^{\mu\nu}`,
:math:`dW` is the normalized distributed random number with mean zero
and variance :math:`\sqrt{\Delta t}`, and :math:`\Delta t` is the time
step for stochastic integration.

Wave-particle interaction
^^^^^^^^^^^^^^^^^^^^^^^^^

For a 2D problem,
reference [Skilling75]_ shows that for forward and
backward propagating Alfvén waves,

.. math::

   \begin{aligned}
     \boldsymbol{u} & = \boldsymbol{v}_0 + \left<\frac{3}{2}(1-\mu^2)\frac{\nu^+ - \nu^-}{\nu^+ + \nu^-}\right>,
     \text{ the velocity of mean wave frame} \\
     \kappa_\parallel & = v^2\left<\frac{1-\mu^2}{2(\nu^+ + \nu^-)}\right>,
     \text{ parallel spatial diffusion coefficient} \\
     D_{pp} & = 4\gamma^2m^2v_A^2\left<\frac{1-\mu^2}{2}\frac{\nu^+\nu^-}{\nu^+ + \nu^-}\right>,
     \text{ momentum diffusion coefficient}
   \end{aligned}

where :math:`\left<\dots\right>` donates :math:`\mu`-average,
:math:`\nu^+` and :math:`\nu^-` are collision frequency against forward
waves and backward waves, respectively. If :math:`\nu^+` is equal to
:math:`\nu^-`,

.. math:: D_{pp} = \frac{1}{9}\frac{p^2v_A^2}{\kappa_\parallel}

where :math:`p=\gamma mv` is particle momentum. Depending on the plasma
parameter and wave properties, we may have to use more complicated
models [Schlickeiser89]_ [Schlickeiser98]_ [LeRoux07]_. The corresponding SDE is

.. math::

   \begin{aligned}
     dp & = \left(-\frac{p}{3}\nabla\cdot\boldsymbol{V} + \frac{4pv_A^2}{9\kappa_\parallel}\right)ds +
     \sqrt{\frac{2p^2v_A^2}{9\kappa_\parallel}}dW(s), \text{if $\kappa_\parallel$ is independent of $p$}\\
     dp & = \left(-\frac{p}{3}\nabla\cdot\boldsymbol{V} + \frac{8pv_A^2}{27\kappa_\parallel}\right)ds +
     \sqrt{\frac{2p^2v_A^2}{9\kappa_\parallel}}dW(s), \text{if $\kappa_\parallel\sim p^{4/3}$}
   \end{aligned}

which are normalized to

.. math::

   \begin{aligned}
     d\tilde{p}_n & = \left(-\frac{\tilde{p}_n}{3}\tilde{\nabla}\cdot\tilde{\boldsymbol{V}} + \frac{4\tilde{p}_n\tilde{v}_A^2}{9\tilde{\kappa}_\parallel}\right)d\tilde{s} + \tilde{p}_n\tilde{v}_A\sqrt{\frac{2}{9\tilde{\kappa}_\parallel}}dW(\tilde{s}), \text{if $\kappa_\parallel$ is independent of $p$}\\
     d\tilde{p}_n & = \left(-\frac{\tilde{p}_n}{3}\tilde{\nabla}\cdot\tilde{\boldsymbol{V}} + \frac{8\tilde{p}_n\tilde{v}_A^2}{27\tilde{\kappa}_\parallel}\right)d\tilde{s} + \tilde{p}_n\tilde{v}_A\sqrt{\frac{2}{9\tilde{\kappa}_\parallel}}dW(\tilde{s}), \text{if $\kappa_\parallel\sim p^{4/3}$}
   \end{aligned}

where :math:`\tilde{p}_n=\tilde{p}\tilde{p}_{n0}=p\tilde{p}_{n0}/p_0`,
where is :math:`\tilde{p}_{n0}` is the numerical value for particles
with :math:`p_0` in the code (e.g., 0.1 as often used),
:math:`\tilde{\nabla}=L_0\nabla`,
:math:`\tilde{\boldsymbol{V}}=\boldsymbol{V}/v_{A0}`,
:math:`\tilde{v}_A=\tilde{v}_{A0}`,
:math:`\tilde{\kappa}_\parallel=\kappa_\parallel/\kappa_0`,
:math:`\kappa_0=L_0v_{A0}`, :math:`\tilde{s}=s/t_0`, and
:math:`t_0=L_0/v_{A0}`. These are all given in the code.

Flow shear
^^^^^^^^^^

For isotropic particle distributions, the flow shear introduces another
momentum diffusion term. If there is no average magnetic
field [Earl88]_.

.. math::

   \begin{aligned}
     D_{pp} & = \Gamma\tau p^2, \\
     \Gamma & = \frac{1}{30}\left(\frac{\partial U_i}{\partial x_j} +
     \frac{\partial U_j}{\partial x_i}\right)^2 -
     \frac{2}{45}\frac{\partial U_i}{\partial x_i}\frac{\partial U_j}{\partial x_j}
     = \frac{2}{15}\sum_{ij}\sigma_{ij}^2
   \end{aligned}

where :math:`\Gamma` is the coefficient of viscous momentum transfer,
:math:`\sigma_{ij}=(\partial_iU_j + \partial_jU_i - 2\nabla\cdot\boldsymbol{U}\delta_{ij}/3)/2`
is the shear tensor, :math:`\tau` is the relaxation time for particle
scattering. According to [Webb18]_, :math:`\tau`
is related particle diffusion coefficient :math:`\kappa_\parallel=v^2\tau/3`. The corresponding SDE is

.. math::

   \begin{aligned}
     dp = \left(-\frac{p}{3}\nabla\cdot\boldsymbol{v} +
     \frac{\Gamma}{p^2}\frac{\partial(p^4\tau)}{\partial p}\right)ds +
     \sqrt{2\Gamma\tau p^2}dW(s)
   \end{aligned}

For :math:`\tau\sim\tau_0(p_0/p)^\alpha`,

.. math::

   \begin{aligned}
     dp = \left(-\frac{p}{3}\nabla\cdot\boldsymbol{v} + \frac{\Gamma\tau_0p_0^\alpha}{p^2}(4-\alpha)p^{3-\alpha}\right)ds +
     \sqrt{2\Gamma\tau_0 p_0^\alpha p^{2-\alpha}}dW(s)
   \end{aligned}

which is normalized to

.. math::

   \begin{aligned}
     d\tilde{p}_n = \left(-\frac{\tilde{p}_n}{3}\tilde{\nabla}\cdot\tilde{\boldsymbol{v}} + (4-\alpha)\tilde{\Gamma}\tilde{\tau}_0\tilde{p}_n^{1-\alpha}\tilde{p}_{n0}^\alpha\right)d\tilde{s} + \sqrt{2\tilde{\Gamma}\tilde{\tau}_0\tilde{p}_n^{2-\alpha}\tilde{p}_{n0}^\alpha}dW(\tilde{s})
   \end{aligned}

where :math:`\tilde{p}_n=\tilde{p}\tilde{p}_{n0}=p\tilde{p}_{n0}/p_0`,
where is :math:`\tilde{p}_{n0}` is the numerical value for particles
with :math:`p_0` in the code (e.g., 0.1 as often used),
:math:`\tilde{\nabla}=L_0\nabla`,
:math:`\tilde{\boldsymbol{v}}=\boldsymbol{v}/v_{A0}`,
:math:`\tilde{\Gamma}=\Gamma t_0^2`, :math:`\tilde{\tau}_0=\tau_0/t_0`,
:math:`\tilde{s}=s/t_0`, and :math:`t_0=L_0/v_{A0}`. For
:math:`\tau\sim\tau_0(p_0/p)^2` [Earl88]_,

.. math::

   \begin{aligned}
     d\tilde{p}_n = \left(-\frac{\tilde{p}_n}{3}\tilde{\nabla}\cdot\tilde{\boldsymbol{v}} + \frac{2\tilde{\Gamma}\tilde{\tau}_0\tilde{p}_{n0}^2}{\tilde{p}_n}\right)d\tilde{s} + \sqrt{2\tilde{\Gamma}\tilde{\tau}_0\tilde{p}_{n0}^2}dW(\tilde{s})
   \end{aligned}

For
:math:`\tau\sim\tau_0(p_0/p)^{2/3}` [Giacalone99]_,

.. math::

   \begin{aligned}
     d\tilde{p}_n & = \left(-\frac{\tilde{p}_n}{3}\tilde{\nabla}\cdot\tilde{\boldsymbol{v}} + \frac{10}{3}\tilde{\Gamma}\tilde{\tau}_0\tilde{p}_{n}^{1/3}\tilde{p}_{n0}^{2/3}\right)d\tilde{s} + \sqrt{2\tilde{\Gamma}\tilde{\tau}_0\tilde{p}_n^{4/3}\tilde{p}_{n0}^{2/3}}dW(\tilde{s}) \\
     \tau_0 & = 3\kappa_{\parallel 0} / v_0^2
   \end{aligned}

If there is an average magnetic field, the equation is more complicated
(see [Williams91]_ [Williams93]_).

.. [Earl88] Earl, J.A., Jokipii, J.R. and Morfill, G., 1988. Cosmic-ray viscosity. The Astrophysical Journal, 331, pp.L91-L94.
.. [LeRoux07] Le Roux, J.A. and Webb, G.M., 2007. Nonlinear cosmic-ray diffusive transport in combined two-dimensional and slab magnetohydrodynamic turbulence: a BGK-Boltzmann approach. The Astrophysical Journal, 667(2), p.930.
.. [Schlickeiser89] Schlickeiser, R., 1989. Cosmic-ray transport and acceleration. I-Derivation of the kinetic equation and application to cosmic rays in static cold media. II-Cosmic rays in moving cold media with application to diffusive shock wave acceleration. The Astrophysical Journal, 336, pp.243-293.
.. [Schlickeiser98] Schlickeiser, R. and Miller, J.A., 1998. Quasi-linear theory of cosmic ray transport and acceleration: the role of oblique magnetohydrodynamic waves and transit-time damping. The Astrophysical Journal, 492(1), p.352.
.. [Skilling75] Skilling, J., 1975. Cosmic Ray Streaming—II effect of particles on alfvén waves. Monthly Notices of the Royal Astronomical Society, 173(2), pp.245-254.
.. [Webb18] Webb, G. M., Barghouty, A. F., Hu, Q., & le Roux, J. A. 2018, The Astrophysical Journal, 855, 31
.. [Williams91] Williams, L.L. and Jokipii, J.R., 1991. Viscosity and inertia in cosmic-ray transport-Effects of an average magnetic field. The Astrophysical Journal, 371, pp.639-647.
.. [Williams93] Williams, L.L., Schwadron, N., Jokipii, J.R. and Gombosi, T.I., 1993. A unified transport equation for both cosmic rays and thermal particles. The Astrophysical Journal, 405, pp.L79-L81.

1D & 2D Models
--------------

In a 1D model, the Fokker-Planck form of the transport is equivalent to the SDEs

.. math::

   dX & = \left(\frac{\partial\kappa}{\partial x} + V_x\right)ds +
     \sqrt{2\kappa} dW_\sigma(s),\\
   dp & = -\frac{p}{3}\frac{\partial V_x}{\partial x}ds,

.. note:: 
  In a 1D or 2D model, the particle drift :math:`\boldsymbol{V}_d` is out of the simulation plane and is not considered here. 

where :math:`dW_\sigma` is the normalized distributed random number with mean
zero and variance :math:`\sqrt{\Delta t}`, and :math:`\Delta t` is the
time step for stochastic integration. This corresponds to a Wiener
process. Numerical approximation is often used for the Wiener process to
replace the normal distribution. We use a uniform distribution in
:math:`[-\sqrt{3}, \sqrt{3}]` in the code.

.. note:: 
  
  Future tests should be included to test this method thoroughly.

In a 2D model, the corresponding SDEs are

.. math::

   dX & = (\nabla\cdot\boldsymbol{\kappa} + \boldsymbol{V})ds +
     \sum_\sigma\boldsymbol{\alpha}_\sigma dW_\sigma(s),\\
   dp & = -\frac{p}{3}(\nabla\cdot\boldsymbol{V})ds,

where
:math:`\sum_\sigma\alpha_\sigma^\mu\alpha_\sigma^\nu = 2\kappa^{\mu\nu}`.

.. math::

   \boldsymbol{\alpha}_1 =
     \begin{pmatrix}
       \sqrt{2\kappa_\perp} \\
       0
     \end{pmatrix}, \quad
     \boldsymbol{\alpha}_2 =
     \begin{pmatrix}
       0 \\
       \sqrt{2\kappa_\perp}
     \end{pmatrix}, \quad
     \boldsymbol{\alpha}_3 =
     \sqrt{2(\kappa_\parallel - \kappa_\perp)}
     \begin{pmatrix}
       B_x/B \\
       B_y/B
     \end{pmatrix}.

The parameters used at particle locations are calculated from
:math:`v_x`, :math:`v_y`, :math:`B_x`, :math:`B_y`,
:math:`\nabla\cdot\boldsymbol{v}`, :math:`\partial B_x/\partial x`,
:math:`\partial B_x/\partial y`, :math:`\partial B_y/\partial x`, and
:math:`\partial B_y/\partial y`, which are all obtained from the MHD
simulations. We interpolate these parameters to the particle positions
and then calculate the other required parameters:

.. math::

   \begin{aligned}
     \frac{\partial\kappa_{xx}}{\partial x} & = \frac{\partial\kappa_\perp}{\partial x} -
     \frac{\partial(\kappa_\perp-\kappa_\parallel)}{\partial x}\frac{B_x^2}{B^2} -
     2(\kappa_\perp-\kappa_\parallel)\frac{\frac{\partial B_x}{\partial x}B_xB-
     \frac{\partial B}{\partial x}B_x^2}{B^3}, \\
     \frac{\partial\kappa_{yy}}{\partial y} & = \frac{\partial\kappa_\perp}{\partial y} -
     \frac{\partial(\kappa_\perp-\kappa_\parallel)}{\partial y}\frac{B_y^2}{B^2} -
     2(\kappa_\perp-\kappa_\parallel)\frac{\frac{\partial B_y}{\partial y}B_yB-
     \frac{\partial B}{\partial y}B_y^2}{B^3}, \\
     \frac{\partial\kappa_{xy}}{\partial x} & =
     -\frac{\partial(\kappa_\perp-\kappa_\parallel)}{\partial x}
     \frac{B_xB_y}{B^2} - (\kappa_\perp-\kappa_\parallel)
     \frac{\left(\frac{\partial B_x}{\partial x}B_y+
     B_x\frac{\partial B_y}{\partial x}\right)B -
     2B_xB_y\frac{\partial B}{\partial x}}{B^3}, \\
     \frac{\partial\kappa_{xy}}{\partial y} & =
     -\frac{\partial(\kappa_\perp-\kappa_\parallel)}{\partial y}
     \frac{B_xB_y}{B^2} - (\kappa_\perp-\kappa_\parallel)
     \frac{\left(\frac{\partial B_x}{\partial y}B_y+
     B_x\frac{\partial B_y}{\partial y}\right)B -
     2B_xB_y\frac{\partial B}{\partial y}}{B^3}, \\
     \frac{\partial B}{\partial x} & = \frac{1}{B}\left(B_x
     \frac{\partial B_x}{\partial x} + B_y\frac{\partial B_y}{\partial x}\right), \\
     \frac{\partial B}{\partial y} & =
     \frac{1}{B}\left(B_x\frac{\partial B_x}{\partial y} +
     B_y\frac{\partial B_y}{\partial y}\right).
   \end{aligned}

where :math:`\kappa_\parallel` and :math:`\kappa_\perp` can be functions
of :math:`B_x`, :math:`B_y` and :math:`B`, so
:math:`\partial \kappa_\parallel/\partial x`,
:math:`\partial \kappa_\parallel/\partial y`,
:math:`\partial \kappa_\perp/\partial x`, and
:math:`\partial \kappa_\perp/\partial y` still depend on the derivatives
:math:`\partial B_x/\partial x`, :math:`\partial B_x/\partial y`,
:math:`\partial B_y/\partial x`, and :math:`\partial B_y/\partial y`.
The detailed expressions depend on the diffusion model to choose. Considering the expression of :math:`\kappa_\parallel`, we get

.. math::

   \begin{aligned}
     \frac{\partial\kappa}{\partial x}\sim\kappa\left(
     \frac{2}{3L_c}\frac{\partial L_c}{\partial x} -
     \frac{1}{3B}\frac{\partial B}{\partial x} -
     \frac{1}{\sigma^2}\frac{\partial(\sigma^2)}{\partial x}
     \right)
   \end{aligned}

.. note:: 
  We typically assume :math:`L_c` and :math:`\sigma^2` are spatially independent. However, they could vary spatially in a large system.

Time step criteria
^^^^^^^^^^^^^^^^^^

For a 1D problem, the particle moves a distance satisfying :math:`l_x^2=\text{max}\left(\left<\Delta x\right>^2, \left<\Delta x^2\right>\right)` [Strauss17]_, where

.. math::

   \begin{aligned}
     \left<\Delta x\right> = \left(V_x + \frac{d\kappa(x)}{dx}\right)\Delta t,
     \quad \left<\Delta x^2\right> = 2\kappa(x)\Delta t,
   \end{aligned}

and :math:`l_x` should be much smaller than the spatial variation scale
of the fields. In this code, we assume
:math:`\left<\Delta x\right>^2 < \left<\Delta x^2\right>` and choose
:math:`\Delta t` so that :math:`l_x\ll\delta_x`, where :math:`\delta_x`
is the grid size. For the 2D problems, we choose the following criteria
to determine the time step:

.. math::

   \begin{aligned}
     \Delta t_x & = \text{min}\left[\frac{\delta x}{|V_x + \partial_x\kappa_{xx} + \partial_y\kappa_{xy}|},
     \frac{2\kappa_\perp} {(V_x + \partial_x\kappa_{xx} + \partial_y\kappa_{xy})^2}\right], \\
     \Delta t_y & = \text{min}\left[\frac{\delta y}{|V_y + \partial_y\kappa_{yy} +
     \partial_x\kappa_{xy}|},
     \frac{2\kappa_\perp}{(V_y + \partial_y\kappa_{yy} + \partial_x\kappa_{xy})^2}\right],\\
     \Delta t & = \text{min}(\Delta t_x, \Delta t_y).
   \end{aligned}

Higher-order method
^^^^^^^^^^^^^^^^^^^

To get a higher-order solution, we can use a derivative-free Milstein method [Burrage04]_ to solve the SDEs. It is different from the usual method due to one more term, which makes it a higher-order method.

.. note:: 

  This method was used by [Li18]_ but has since not been supported in default.

.. math::

   \begin{aligned}
     dX_t & = f(X_t,t)dt + g(X_t,t)dW_t, \\
     X_{n+1} & = X_n + f_n h + g_n\Delta W_n +
     \frac{1}{2\sqrt{h}}[g(\bar{X}_n)-g_n][(\Delta W_n)^2-h], \\
     \bar{X}_n & = X_n + f_n h + g_n\sqrt{h}, \\
     \Delta W_n & = [W_{t+h}-W_t] \sim \sqrt{h}N(0,1),
   \end{aligned}

where :math:`X` corresponds to spatial positions :math:`x`, :math:`y`
and particle momentum :math:`p` in our simulation. Here :math:`f(X_t,t)`
is the deterministic term, :math:`g(X_t,t)` is the probabilistic term,
:math:`h` is the time step, and :math:`N(0,1)` indicates a normal
distribution, which is substituted with a uniform distribution
:math:`[-\sqrt{3}, \sqrt{3}]` in our simulations to speed up the
computation.

.. [Burrage04] Burrage, K., Burrage, P.M. and Tian, T., 2004. Numerical methods for strong solutions of stochastic differential equations: an overview. Proceedings of the Royal Society of London. Series A: Mathematical, Physical and Engineering Sciences, 460(2041), pp.373-402.
.. [Li18] Large-scale Compression Acceleration during Magnetic Reconnection in a Low-β Plasma, Xiaocan Li, Fan Guo, Hui Li, and Shengtai Li, The Astrophysical Journal 866, no. 1 (2018): 4.
.. [Strauss17] Strauss, R. and Effenberger, F., 2017. A hitch-hiker’s guide to stochastic differential equations. Space Science Reviews, 212(1), pp.151-192.

Cartesian Coordinates
---------------------

The Fokker-Planck form of the focused transport equation is

.. math::

   \begin{aligned}
     \frac{\partial F}{\partial t} = & \nabla\cdot\left[\nabla\cdot(\boldsymbol{\kappa}_\perp\nabla F)\right] - \nabla\cdot\left[(v\mu\boldsymbol{b} + \boldsymbol{V} + \boldsymbol{V}_d + \nabla\cdot\boldsymbol{\kappa}_\perp)F\right] \nonumber\\
     & + \frac{\partial}{\partial\mu^2}\left(D_{\mu\mu}F\right) - \frac{\partial}{\partial\mu}\left[\left(\frac{d\mu}{dt}+\frac{\partial D_{\mu\mu}}{\partial\mu}\right)F\right] - \frac{\partial}{\partial p}\left(\frac{dp}{dt}F\right)
   \end{aligned}

where :math:`F=fp^2`. The corresponding SDEs are

.. math::

   \begin{aligned}
     d\boldsymbol{X} & = (v\mu\boldsymbol{b} + \boldsymbol{V} + \boldsymbol{V}_d + \nabla\cdot\boldsymbol{\kappa}_\perp)dt + \sum_\sigma\boldsymbol{\alpha}_\sigma dW_\sigma(s) \\
     dp & = \left(\frac{dp}{dt}\right)dt \\
     d\mu & = \left(\frac{d\mu}{dt}+\frac{\partial D_{\mu\mu}}{\partial\mu}\right)dt + \sqrt{2D_{\mu\mu}}dW_\mu(t)
   \end{aligned}

where
:math:`\sum_\sigma\alpha_\sigma^\mu\alpha_\sigma^\nu = 2\kappa_\perp^{\mu\nu}`.
We can use the results from the 3D model with whole spatial
diffusion tensor :math:`\boldsymbol{\kappa}` but set
:math:`\kappa_\parallel=0`. The matrix

.. math::

   P =
     \begin{pmatrix}
       0 & -b_xb_z\sqrt{2\kappa_\perp}/\sqrt{b_x^2+b_y^2} &
       -b_y\sqrt{2\kappa_\perp}/\sqrt{b_x^2+b_y^2}\\
       0 & -b_yb_z\sqrt{2\kappa_\perp}/\sqrt{b_x^2+b_y^2} &
       b_x\sqrt{2\kappa_\perp}/\sqrt{b_x^2+b_y^2}\\
       0 & \sqrt{b_x^2+b_y^2}\sqrt{2\kappa_\perp}         & 0
     \end{pmatrix}

Since the first column is all zeros, we only need two Weiner processes
to describe the spatial diffusion. For 2D problems,

.. math::

   P = \frac{\sqrt{2\kappa_\perp}}{\sqrt{b_x^2+b_y^2}}
     \begin{pmatrix}
       -b_xb_z &
       -b_y\\
       -b_yb_z &
       b_x
     \end{pmatrix}

To calculate the drift velocity

.. math::

   \begin{aligned}
     & \boldsymbol{V}_d=\frac{cpv}{qB}\left\{\frac{1-\mu^2}{2}\frac{\boldsymbol{B}\times\nabla B}{B^2}+\mu^2\frac{\boldsymbol{B}\times[(\boldsymbol{B}\cdot\nabla)\boldsymbol{B}]}{B^3}+\frac{1-\mu^2}{2}\frac{\boldsymbol{B}(\boldsymbol{B}\cdot\nabla\times\boldsymbol{B})}{B^3}\right\},
   \end{aligned}

we need

.. math::

   \begin{aligned}
     (\boldsymbol{B}\times\nabla B)_x & = B_y\partial_z B - B_z\partial_y B \\
     (\boldsymbol{B}\times\nabla B)_y & = B_z\partial_x B - B_x\partial_z B \\
     (\boldsymbol{B}\times\nabla B)_z & = B_x\partial_y B - B_y\partial_x B \\
     \{\boldsymbol{B}\times[(\boldsymbol{B}\cdot\nabla)\boldsymbol{B}]\}_x & = B_y(\boldsymbol{B}\cdot\nabla)B_z - B_z(\boldsymbol{B}\cdot\nabla)B_y \\
     \{\boldsymbol{B}\times[(\boldsymbol{B}\cdot\nabla)\boldsymbol{B}]\}_y & = B_z(\boldsymbol{B}\cdot\nabla)B_x - B_x(\boldsymbol{B}\cdot\nabla)B_z \\
     \{\boldsymbol{B}\times[(\boldsymbol{B}\cdot\nabla)\boldsymbol{B}]\}_z & = B_x(\boldsymbol{B}\cdot\nabla)B_y - B_y(\boldsymbol{B}\cdot\nabla)B_x \\
     \boldsymbol{B}\cdot\nabla & = B_x\partial_x + B_y\partial_y + B_z\partial_z \\
     \boldsymbol{B}\cdot\nabla\times\boldsymbol{B} & = B_x (\partial_y B_z - \partial_z B_y) + \nonumber \\
     & B_y (\partial_z B_x - \partial_x B_z) + B_z (\partial_x B_y - \partial_y B_x)
   \end{aligned}

For 2D problems, these can be simplified using

.. math::

   \begin{aligned}
     (\boldsymbol{B}\times\nabla B)_x & = - B_z\partial_y B \\
     (\boldsymbol{B}\times\nabla B)_y & = B_z\partial_x B \\
     (\boldsymbol{B}\times\nabla B)_z & = B_x\partial_y B - B_y\partial_x B \\
     \boldsymbol{B}\cdot\nabla & = B_x\partial_x + B_y\partial_y \\
     \boldsymbol{B}\cdot\nabla\times\boldsymbol{B} & = B_x\partial_y B_z - B_y\partial_x B_z + B_z (\partial_x B_y - \partial_y B_x)
   \end{aligned}

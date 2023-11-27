Spherical Coordinates
---------------------

The spatial diffusion coefficient is in the same form.

.. math::

   \begin{aligned}
     & \boldsymbol{\kappa}_\perp =
     \begin{bmatrix}
       \kappa_{\perp rr} & \kappa_{\perp r\theta} & \kappa_{\perp r\phi} \\
       \kappa_{\perp r\theta} & \kappa_{\perp\theta\theta} & \kappa_{\perp\theta\phi} \\
       \kappa_{\perp r\phi} & \kappa_{\perp\theta\phi} & \kappa_{\perp\phi\phi}
     \end{bmatrix}
   \end{aligned}

where
:math:`\kappa_{\perp ij}=\kappa_\perp\delta_{ij} - \kappa_\perp b_ib_j`,
and :math:`i,j` are :math:`r,\theta,\phi`. The gradients of
:math:`\kappa_{ij}` are

.. math::

   \begin{aligned}
     \partial_i\kappa_{\perp ij} & = \delta_{ij}\partial_i\kappa_\perp -b_ib_j\partial_i\kappa_\perp -
     \frac{\kappa_\perp}{B}(b_j\partial_iB_i + b_i\partial_iB_j -
     2b_ib_j\partial_iB)
   \end{aligned}

We then need to transfer the focused transport equation to the spherical
coordinates. Since we don’t have cross-diffusion terms (spatial and
momentum), we can ignore the momentum diffusion for now.

.. math::

   \begin{aligned}
     \frac{\partial F}{\partial t} = & \nabla\cdot(\boldsymbol{\kappa}_\perp\cdot\nabla F) - \nabla\cdot\left[(v\mu\boldsymbol{b} + \boldsymbol{V} + \boldsymbol{V}_d)F\right] \nonumber\\
     & + \frac{\partial}{\partial\mu^2}\left(D_{\mu\mu}F\right) - \frac{\partial}{\partial\mu}\left[\left(\frac{d\mu}{dt}+\frac{\partial D_{\mu\mu}}{\partial\mu}\right)F\right] - \frac{\partial}{\partial p}\left(\frac{dp}{dt}F\right)
   \end{aligned}

The rest is essentially the same as before except for
:math:`\boldsymbol{\kappa}_\perp`, :math:`v\mu\boldsymbol{b}`, and the
pitch-angle diffusion terms. We will first need to change :math:`F` to
:math:`F_1=F\sin\theta r^2` [Jokipii77]_ [Pei10]_.
The resulting SDEs are

.. math::

   \begin{aligned}
     dr & = \frac{dr}{dt}dt + [P.dW_t]_r \\
     d\theta & = \frac{d\theta}{dt}dt + [P.dW_t]_\theta \\
     d\phi & = \frac{d\phi}{dt}dt + [P.dW_t]_\phi\\
     dp & = \frac{dp}{dt}dt \\
     d\mu & = \left(\frac{d\mu}{dt}+\frac{\partial D_{\mu\mu}}{\partial\mu}\right)dt + \sqrt{2D_{\mu\mu}}dW_\mu(t)
   \end{aligned}

where

.. math::

   \begin{aligned}
     \frac{dr}{dt} & = v\mu b_r + V_r+V_{dr} +
     \frac{\partial\kappa_{\perp rr}}{\partial r} + \frac{2}{r}\kappa_{\perp rr}+
     \frac{1}{r}\frac{\partial\kappa_{\perp r\theta}}{\partial\theta}+
     \frac{\cos\theta}{r\sin\theta}\kappa_{\perp r\theta}+
     \frac{1}{r\sin\theta}\frac{\partial\kappa_{\perp r\phi}}{\partial\phi} \\
     \frac{d\theta}{dt} & = \frac{v\mu b_\theta + V_\theta+V_{d\theta}}{r} +
     \frac{1}{r}\frac{\partial\kappa_{\perp r\theta}}{\partial r} + \frac{\kappa_{\perp r\theta}}{r^2}+
     \frac{1}{r^2}\frac{\partial\kappa_{\perp \theta\theta}}{\partial\theta}+
     \frac{\cos\theta}{r^2\sin\theta}\kappa_{\perp \theta\theta}+
     \frac{1}{r^2\sin\theta}\frac{\partial\kappa_{\perp \theta\phi}}{\partial\phi} \\
     \frac{d\phi}{dt} & = \frac{v\mu b_\phi + V_\phi+V_{d\phi}}{r\sin\theta}+
     \frac{1}{r\sin\theta}\frac{\partial\kappa_{\perp r\phi}}{\partial r} +
     \frac{\kappa_{\perp r\phi}}{r^2\sin\theta} +
     \frac{1}{r^2\sin\theta}\frac{\partial\kappa_{\perp \theta\phi}}{\partial\theta}+
     \frac{1}{r^2\sin^2\theta}\frac{\partial\kappa_{\perp \phi\phi}}{\partial\phi}
   \end{aligned}

To calculate :math:`d\mu/dt` and :math:`dp/dt`, we will need

.. math::

   \begin{aligned}
     -\boldsymbol{b}\cdot\nabla\ln B & = -\frac{1}{B}\left(b_r\partial_r B + \frac{b_\theta}{r}\partial_\theta B + \frac{b_\phi}{r\sin\theta}\partial_\phi B\right) \\
     \nabla\cdot\boldsymbol{V} & = \partial_r V_r +\frac{2V_r}{r}+
     \frac{1}{r}\partial_\theta V_\theta +
     \frac{\cos\theta}{r\sin\theta}V_\theta +
     \frac{1}{r\sin\theta}\partial_\phi V_\phi
   \end{aligned}

:math:`b_ib_j\frac{\partial V_i}{\partial x_j}` is more complicated.

.. math::

   \begin{aligned}
     b_ib_j\frac{\partial V_i}{\partial x_j} =
     & b_r^2\partial_r V_r + \frac{b_rb_\theta}{r}\partial_\theta V_r + \frac{b_rb_\phi}{r\sin\theta}\partial_\phi V_r - \frac{b_rb_\theta V_\theta + b_rb_\phi V_\phi}{r} + \nonumber\\
     & b_rb_\theta\partial_r V_\theta + \frac{b_\theta^2}{r}\partial_\theta V_\theta + \frac{b_\theta b_\phi}{r\sin\theta}\partial_\phi V_\theta + \frac{b_\theta^2V_r}{r} - \frac{\cot\theta b_\theta b_\phi V_\phi}{r} + \nonumber\\
     & b_rb_\phi\partial_r V_\phi + \frac{b_\theta b_\phi}{r}\partial_\theta V_\phi + \frac{b_\phi^2}{r\sin\theta}\partial_\phi V_\phi + \frac{b_\phi^2V_r}{r} + \frac{\cot\theta b_\phi^2V_\theta}{r}
   \end{aligned}

Similar for :math:`b_iV_j\frac{\partial V_i}{\partial x_j}`,

.. math::

   \begin{aligned}
     b_iV_j\frac{\partial V_i}{\partial x_j} =
     & b_rV_r\partial_r V_r + \frac{b_rV_\theta}{r}\partial_\theta V_r + \frac{b_rV_\phi}{r\sin\theta}\partial_\phi V_r - \frac{b_r(V_\theta^2+V_\phi^2)}{r} + \nonumber\\
     & b_\theta V_r\partial_r V_\theta + \frac{b_\theta V_\theta}{r}\partial_\theta V_\theta + \frac{b_\theta V_\phi}{r\sin\theta}\partial_\phi V_\theta + \frac{b_\theta V_\theta V_r}{r} - \frac{\cot\theta b_\theta V_\phi^2}{r} + \nonumber\\
     & b_\phi V_r\partial_r V_\phi + \frac{b_\phi V_\theta}{r}\partial_\theta V_\phi + \frac{b_\phi V_\phi}{r\sin\theta}\partial_\phi V_\phi + \frac{b_\phi V_\phi V_r}{r} + \frac{\cot\theta b_\phi V_\phi V_\theta}{r}
   \end{aligned}

For 1D probelms, :math:`F_1=fp^2r^2`, and the corresponding SDE for
:math:`r` is

.. math::

   \begin{aligned}
     dr = & \left(v\mu b_r + V_r +
     \frac{\partial\kappa_{\perp rr}}{\partial r} + \frac{2}{r}\kappa_{\perp rr}\right)dt + \sqrt{2\kappa_{\perp rr}}dW_t
   \end{aligned}

For :math:`d\mu` and :math:`dp`, we will need

.. math::

   \begin{aligned}
     -\boldsymbol{b}\cdot\nabla\ln B = & -\frac{1}{B}b_r\partial_r B \\
     \nabla\cdot\boldsymbol{V} = & \partial_r V_r +\frac{2V_r}{r} \\
     b_ib_j\frac{\partial V_i}{\partial x_j} =
     & b_r^2\partial_r V_r - \frac{b_rb_\theta V_\theta + b_rb_\phi V_\phi}{r} + \nonumber \\
     & b_rb_\theta\partial_r V_\theta + \frac{b_\theta^2V_r}{r} - \frac{\cot\theta b_\theta b_\phi V_\phi}{r} + \nonumber \\
     & b_rb_\phi\partial_r V_\phi + \frac{b_\phi^2V_r}{r} + \frac{\cot\theta b_\phi^2V_\theta}{r} \\
     b_iV_j\frac{\partial V_i}{\partial x_j} =
     & b_rV_r\partial_r V_r - \frac{b_r(V_\theta^2+V_\phi^2)}{r} + \nonumber \\
     & b_\theta V_r\partial_r V_\theta + \frac{b_\theta V_\theta V_r}{r} - \frac{\cot\theta b_\theta V_\phi^2}{r} + \nonumber \\
     & b_\phi V_r\partial_r V_\phi + \frac{b_\phi V_\phi V_r}{r} + \frac{\cot\theta b_\phi V_\phi V_\theta}{r}
   \end{aligned}

For 2D problems (:math:`r-\theta` plane), :math:`F_1=fp^2r^2\sin\theta`,
and the corresponding SDEs are

.. math::

   \begin{aligned}
     dr & = \left(v\mu b_r + V_r+
     \frac{\partial\kappa_{\perp rr}}{\partial r} + \frac{2}{r}\kappa_{\perp rr}+
     \frac{1}{r}\frac{\partial\kappa_{\perp r\theta}}{\partial\theta}+
     \frac{\cos\theta}{r\sin\theta}\kappa_{\perp r\theta}\right)dt + [P.dW_t]_r \\
     d\theta & = \left(\frac{v\mu b_\theta + V_\theta}{r} +
     \frac{1}{r}\frac{\partial\kappa_{\perp r\theta}}{\partial r} + \frac{\kappa_{\perp r\theta}}{r^2}+
     \frac{1}{r^2}\frac{\partial\kappa_{\perp \theta\theta}}{\partial\theta}+
     \frac{\cos\theta}{r^2\sin\theta}\kappa_{\perp \theta\theta}\right)dt + [P.dW_t]_\theta
   \end{aligned}

where :math:`P` has the same form as the one in Parker transport. For
:math:`d\mu` and :math:`dp`, we will need

.. math::

   \begin{aligned}
     -\boldsymbol{b}\cdot\nabla\ln B = & -\frac{1}{B}\left(b_r\partial_r B + \frac{b_\theta}{r}\partial_\theta B\right) \\
     \nabla\cdot\boldsymbol{V} = & \partial_r V_r +\frac{2V_r}{r}+
     \frac{1}{r}\partial_\theta V_\theta +
     \frac{\cos\theta}{r\sin\theta}V_\theta \\
     b_ib_j\frac{\partial V_i}{\partial x_j} =
     & b_r^2\partial_r V_r + \frac{b_rb_\theta}{r}\partial_\theta V_r - \frac{b_rb_\theta V_\theta + b_rb_\phi V_\phi}{r} + \nonumber\\
     & b_rb_\theta\partial_r V_\theta + \frac{b_\theta^2}{r}\partial_\theta V_\theta + \frac{b_\theta^2V_r}{r} - \frac{\cot\theta b_\theta b_\phi V_\phi}{r} + \nonumber\\
     & b_rb_\phi\partial_r V_\phi + \frac{b_\theta b_\phi}{r}\partial_\theta V_\phi + \frac{b_\phi^2V_r}{r} + \frac{\cot\theta b_\phi^2V_\theta}{r} \\
     b_iV_j\frac{\partial V_i}{\partial x_j} =
     & b_rV_r\partial_r V_r + \frac{b_rV_\theta}{r}\partial_\theta V_r - \frac{b_r(V_\theta^2+V_\phi^2)}{r} + \nonumber\\
     & b_\theta V_r\partial_r V_\theta + \frac{b_\theta V_\theta}{r}\partial_\theta V_\theta + \frac{b_\theta V_\theta V_r}{r} - \frac{\cot\theta b_\theta V_\phi^2}{r} + \nonumber\\
     & b_\phi V_r\partial_r V_\phi + \frac{b_\phi V_\theta}{r}\partial_\theta V_\phi + \frac{b_\phi V_\phi V_r}{r} + \frac{\cot\theta b_\phi V_\phi V_\theta}{r}
   \end{aligned}

The drift velocity

.. math::

   \begin{aligned}
     & \boldsymbol{V}_d=\frac{cpv}{qB}\left\{\frac{1-\mu^2}{2}\frac{\boldsymbol{B}\times\nabla B}{B^2}+\mu^2\frac{\boldsymbol{B}\times[(\boldsymbol{B}\cdot\nabla)\boldsymbol{B}]}{B^3}+\frac{1-\mu^2}{2}\frac{\boldsymbol{B}(\boldsymbol{B}\cdot\nabla\times\boldsymbol{B})}{B^3}\right\}
   \end{aligned}

Calculations needed for the first term (gradient drift):

.. math::

   \begin{aligned}
     & (\nabla B)_r=\partial_r B;\quad
     (\nabla B)_\theta=\frac{\partial_\theta B}{r};\quad
     (\nabla B)_\phi=\frac{\partial_\phi B}{r\sin\theta} \\
     & (\boldsymbol{B}\times\nabla B)_r = B_\theta(\nabla B)_\phi - B_\phi(\nabla B)_\theta \\
     & (\boldsymbol{B}\times\nabla B)_\theta = B_\phi(\nabla B)_r - B_r(\nabla B)_\phi \\
     & (\boldsymbol{B}\times\nabla B)_\phi = B_r(\nabla B)_\theta - B_\theta(\nabla B)_r
   \end{aligned}

Calculations needed for the second term (curvature drift):

.. math::

   \begin{aligned}
     & \boldsymbol{C} = (\boldsymbol{B}\cdot\nabla)\boldsymbol{B} \\
     & C_r = B_r\partial_r B_r + \frac{B_\theta}{r}\partial_\theta B_r + \frac{B_\phi}{r\sin\theta}\partial_\phi B_r - \frac{B_\theta^2+B_\phi^2}{r} \\
     & C_\theta = B_r\partial_r B_\theta + \frac{B_\theta}{r}\partial_\theta B_\theta + \frac{B_\phi}{r\sin\theta}\partial_\phi B_\theta + \frac{B_r B_\theta}{r} - \frac{\cot\theta B_\phi^2}{r} \\
     & C_\phi = B_r\partial_r B_\phi + \frac{B_\theta}{r}\partial_\theta B_\phi + \frac{B_\phi}{r\sin\theta}\partial_\phi B_\phi + \frac{B_r B_\phi}{r} + \frac{\cot\theta B_\theta B_\phi}{r} \\
     & (\boldsymbol{B}\times\boldsymbol{C})_r = B_\theta C_\phi - B_\phi C_\theta \\
     & (\boldsymbol{B}\times\boldsymbol{C})_\theta = B_\phi C_r - B_r C_\phi \\
     & (\boldsymbol{B}\times\boldsymbol{C})_\phi = B_r C_\theta - B_\theta C_r
   \end{aligned}

Calculations needed for the third term (parallel drift):

.. math::

   \begin{aligned}
     & (\nabla\times\boldsymbol{B})_r =
     \frac{\partial_\theta(\sin\theta B_\phi)}{r\sin\theta}
     - \frac{\partial_\phi B_\theta}{r\sin\theta} =
     \frac{\partial_\theta B_\phi}{r} +
     \frac{\cos\theta}{r\sin\theta}B_\phi
     - \frac{\partial_\phi B_\theta}{r\sin\theta} \\
     & (\nabla\times\boldsymbol{B})_\theta =
     \frac{\partial_\phi B_r}{r\sin\theta}
     - \frac{\partial_r(rB_\phi)}{r} =
     \frac{\partial_\phi B_r}{r\sin\theta}
     -\partial_r B_\phi - \frac{B_\phi}{r} \\
     & (\nabla\times\boldsymbol{B})_\phi =
     \frac{\partial_r(rB_\theta)}{r}
     - \frac{\partial_\theta B_r}{r} =
     \partial_r B_\theta + \frac{B_\theta}{r}
     - \frac{\partial_\theta B_r}{r}
   \end{aligned}

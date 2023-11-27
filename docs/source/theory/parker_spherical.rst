Spherical Coordinates
---------------------

In spherical coordinates, the drift velocity

.. math::

   \begin{aligned}
     & \boldsymbol{V}_d = \frac{pcw}{3q}\nabla\times\left(\frac{\boldsymbol{B}}{B^2}\right)
     = \frac{1}{3q}\frac{p^2c^3}{\sqrt{p^2c^2+m^2c^4}}
     \left(\frac{1}{B^2}\nabla\times\boldsymbol{B} -
     \frac{2}{B^3}\nabla B\times\boldsymbol{B}\right) \\
     & (\nabla\times\boldsymbol{B})_r =
     \frac{1}{r\sin\theta}\frac{\partial}{\partial\theta}(\sin\theta B_\phi)
     - \frac{1}{r\sin\theta}\frac{\partial B_\theta}{\partial\phi} =
     \frac{1}{r}\frac{\partial B_\phi}{\partial\theta} +
     \frac{\cos\theta}{r\sin\theta}B_\phi
     - \frac{1}{r\sin\theta}\frac{\partial B_\theta}{\partial\phi} \\
     & (\nabla\times\boldsymbol{B})_\theta =
     \frac{1}{r\sin\theta}\frac{\partial B_r}{\partial\phi}
     - \frac{1}{r}\frac{\partial}{\partial r}(rB_\phi) =
     \frac{1}{r\sin\theta}\frac{\partial B_r}{\partial\phi}
     -\frac{\partial B_\phi}{\partial r} - \frac{B_\phi}{r} \\
     & (\nabla\times\boldsymbol{B})_\phi =
     \frac{1}{r}\frac{\partial}{\partial r}(rB_\theta)
     - \frac{1}{r}\frac{\partial B_r}{\partial\theta} =
     \frac{\partial B_\theta}{\partial r} + \frac{B_\theta}{r}
     - \frac{1}{r}\frac{\partial B_r}{\partial\theta} \\
     & (\nabla B)_r=\frac{\partial B}{\partial r};\quad
     (\nabla B)_\theta=\frac{1}{r}\frac{\partial B}{\partial\theta};\quad
     (\nabla B)_\phi=\frac{1}{r\sin\theta}\frac{\partial B}{\partial\phi} \\
     & (\nabla B\times\boldsymbol{B})_r = (\nabla B)_\theta B_\phi - (\nabla B)_\phi B_\theta \\
     & (\nabla B\times\boldsymbol{B})_\theta = (\nabla B)_\phi B_r - (\nabla B)_r B_\phi \\
     & (\nabla B\times\boldsymbol{B})_\phi = (\nabla B)_r B_\theta - (\nabla B)_\theta B_r
   \end{aligned}

The spatial diffusion coefficient is in the same form.

.. math::

   \begin{aligned}
     & \kappa =
     \begin{bmatrix}
       \kappa_{rr} & \kappa_{r\theta} & \kappa_{r\phi} \\
       \kappa_{r\theta} & \kappa_{\theta\theta} & \kappa_{\theta\phi} \\
       \kappa_{r\phi} & \kappa_{\theta\phi} & \kappa_{\phi\phi}
     \end{bmatrix}
     % & \kappa_{rr} = \kappa_\perp - \frac{\kappa_\perp-\kappa_\parallel}{B^2}B_r^2\\
     % & \kappa_{\theta\theta} = \kappa_\perp - \frac{\kappa_\perp-\kappa_\parallel}{B^2}B_\theta^2\\
     % & \kappa_{\phi\phi} = \kappa_\perp - \frac{\kappa_\perp-\kappa_\parallel}{B^2}B_\phi^2\\
     % & \kappa_{r\theta} = - \frac{\kappa_\perp-\kappa_\parallel}{B^2}B_rB_\theta\\
     % & \kappa_{r\phi} = - \frac{\kappa_\perp-\kappa_\parallel}{B^2}B_rB_\phi\\
     % & \kappa_{\theta\phi} = - \frac{\kappa_\perp-\kappa_\parallel}{B^2}B_\theta B_\phi
   \end{aligned}

where
:math:`\kappa_{ij}=\kappa_\perp\delta_{ij} - (\kappa_\perp - \kappa_\parallel)b_ib_j`,
and :math:`i,j` are :math:`r,\theta,\phi`.

Since :math:`\kappa` is of the same form as that in the Cartesian
coordinates, the gradients of :math:`\kappa` are

.. math::

   \begin{aligned}
     \partial_r\kappa_{rr} & = \partial_r\kappa_\perp +
     \partial_r(\kappa_\parallel-\kappa_\perp)b_r^2 +
     2(\kappa_\parallel-\kappa_\perp)(B_rB\partial_rB_r - B_r^2\partial_r B)/B^3, \\
     \partial_\theta\kappa_{\theta\theta} & = \partial_\theta\kappa_\perp +
     \partial_\theta(\kappa_\parallel-\kappa_\perp)b_\theta^2 +
     2(\kappa_\parallel-\kappa_\perp)(B_\theta B\partial_\theta B_\theta - B_\theta^2\partial_\theta B)/B^3, \\
     \partial_\phi\kappa_{\phi\phi} & = \partial_\phi\kappa_\perp +
     \partial_\phi(\kappa_\parallel-\kappa_\perp)b_\phi^2 +
     2(\kappa_\parallel-\kappa_\perp)(B_\phi B\partial_\phi B_\phi - B_\phi^2\partial_\phi B)/B^3, \\
     \partial_r\kappa_{r\theta} & = \partial_r(\kappa_\parallel-\kappa_\perp)b_rb_\theta +
     (\kappa_\parallel-\kappa_\perp)[(B_\theta\partial_rB_r + B_r\partial_rB_\theta)B -
     2B_rB_\theta\partial_rB] / B^3, \\
     \partial_\theta\kappa_{r\theta} & = \partial_\theta(\kappa_\parallel-\kappa_\perp)b_rb_\theta +
     (\kappa_\parallel-\kappa_\perp)[(B_\theta\partial_\theta B_r + B_r\partial_\theta B_\theta)B -
     2B_rB_\theta\partial_\theta B] / B^3, \\
     \partial_r\kappa_{r\phi} & = \partial_r(\kappa_\parallel-\kappa_\perp)b_rb_\phi +
     (\kappa_\parallel-\kappa_\perp)[(B_\phi\partial_rB_r + B_r\partial_rB_\phi)B -
     2B_rB_\phi\partial_rB] / B^3, \\
     \partial_\phi\kappa_{r\phi} & = \partial_\phi(\kappa_\parallel-\kappa_\perp)b_rb_\phi +
     (\kappa_\parallel-\kappa_\perp)[(B_\phi\partial_\phi B_r + B_r\partial_\phi B_\phi)B -
     2B_rB_\phi\partial_\phi B] / B^3, \\
     \partial_\theta\kappa_{\theta\phi} & = \partial_\theta(\kappa_\parallel-\kappa_\perp)b_\theta b_\phi +
     (\kappa_\parallel-\kappa_\perp)[(B_\phi\partial_\theta B_\theta + B_\theta\partial_\theta B_\phi)B -
     2B_\theta B_\phi\partial_\theta B] / B^3, \\
     \partial_\phi\kappa_{\theta\phi} & = \partial_\phi(\kappa_\parallel-\kappa_\perp)b_\theta b_\phi +
     (\kappa_\parallel-\kappa_\perp)[(B_\phi\partial_\phi B_\theta + B_\theta\partial_\phi B_\phi)B -
     2B_\theta B_\phi\partial_\phi B] / B^3.
   \end{aligned}

We then need to transfer the Parker transport equation to the spherical
coordinates. Since we don’t have cross-diffusion terms (spatial and
momentum), we can ignore the momentum diffusion for now.

.. note:: 
  For a more complete equation, the cross-diffusion terms should be included.

.. math::

   \begin{aligned}
     \frac{\partial F}{\partial t}
     & = -(\boldsymbol{V}+\boldsymbol{V}_d)\cdot\nabla F
     - (\nabla\cdot\boldsymbol{V})F
     + \frac{\partial}{\partial p}\left[\frac{p}{3}(\nabla\cdot\boldsymbol{V})F\right]
     + \nabla\cdot(\boldsymbol{\kappa}\cdot\nabla F)
   \end{aligned}

where :math:`F=fp^2`. Since :math:`\nabla\cdot\boldsymbol{V}_d=0`, we
can add one more term :math:`-(\nabla\cdot\boldsymbol{V}_d)F` to the
right. Then,

.. math::

   \begin{aligned}
     \frac{\partial F}{\partial t}
     & = -\nabla\cdot((\boldsymbol{V}+\boldsymbol{V}_d)F)
     + \frac{\partial}{\partial p}\left[\frac{p}{3}(\nabla\cdot\boldsymbol{V})F\right]
     + \nabla\cdot(\boldsymbol{\kappa}\cdot\nabla F)
   \end{aligned}

Taking :math:`\boldsymbol{V}+\boldsymbol{V}_d\rightarrow\boldsymbol{V}`,

.. math::

   \begin{aligned}
     \nabla\cdot(\boldsymbol{V}F) & =
     \frac{1}{r^2}\frac{\partial}{\partial r}(r^2V_rF)
     +\frac{1}{r\sin\theta}\frac{\partial}{\partial\theta}(\sin\theta V_\theta F)
     +\frac{1}{r\sin\theta}\frac{\partial}{\partial\phi}(V_\phi F) \\
     & = \frac{\partial(V_rF)}{\partial r} + \frac{2}{r}V_rF
     +\frac{\partial}{\partial\theta}\left(\frac{V_\theta F}{r}\right)
     +\frac{\cos\theta}{r\sin\theta}V_\theta F
     +\frac{\partial}{\partial\phi}\left(\frac{V_\phi F}{r\sin\theta}\right)
   \end{aligned}

so there is 2 additional terms (2nd and 4th) if we want to write the
equation Fokker–Planck form. It turns out that we need to change
:math:`F` to
:math:`F_1=F\sin\theta r^2` [Jokipii77]_ [Pei10]_.
Multiplying the above equation by :math:`r^2\sin\theta`, we get

.. math::

   \begin{aligned}
     r^2\sin\theta\nabla\cdot(\boldsymbol{V}F) & =
     \frac{\partial(V_rF_1)}{\partial r}
     +\frac{\partial}{\partial\theta}\left(\frac{V_\theta F_1}{r}\right)
     +\frac{\partial}{\partial\phi}\left(\frac{V_\phi F_1}{r\sin\theta}\right)
   \end{aligned}

For the diffusion term,

.. math::

   \begin{aligned}
     \boldsymbol{\kappa}\cdot\nabla F = &
     \left(\kappa_{rr}\frac{\partial F}{\partial r} +
     \kappa_{r\theta}\frac{1}{r}\frac{\partial F}{\partial\theta} +
     \kappa_{r\phi}\frac{1}{r\sin\theta}\frac{\partial F}{\partial\phi}
     \right)\hat{e}_r + \\\nonumber
     & \left(\kappa_{r\theta}\frac{\partial F}{\partial r} +
     \kappa_{\theta\theta}\frac{1}{r}\frac{\partial F}{\partial\theta} +
     \kappa_{\theta\phi}\frac{1}{r\sin\theta}\frac{\partial F}{\partial\phi}
     \right)\hat{e}_\theta + \\\nonumber
     & \left(\kappa_{r\phi}\frac{\partial F}{\partial r} +
     \kappa_{\theta\phi}\frac{1}{r}\frac{\partial F}{\partial\theta} +
     \kappa_{\phi\phi}\frac{1}{r\sin\theta}\frac{\partial F}{\partial\phi}
     \right)\hat{e}_\phi
   \end{aligned}

Taking :math:`\boldsymbol{A}=\boldsymbol{\kappa}\cdot\nabla F`,

.. math::

   \begin{aligned}
     r^2\sin\theta\nabla\cdot\boldsymbol{A} =
     \frac{\partial(r^2\sin\theta A_r)}{\partial r} +
     \frac{\partial(r\sin\theta A_\theta)}{\partial\theta} +
     \frac{\partial(rA_\phi)}{\partial\phi}
   \end{aligned}

The 1st term on the right is expanded to

.. math::

   \begin{aligned}
     & \frac{\partial^2}{\partial r^2}(\kappa_{rr}F_1) +
     \frac{\partial^2}{\partial r\partial\theta}\left(\frac{\kappa_{r\theta}}{r}F_1\right) +
     \frac{\partial^2}{\partial r\partial\phi}\left(\frac{\kappa_{r\phi}}{r\sin\theta}F_1\right) \\
     \nonumber
     & -\frac{\partial}{\partial r}\left[\left(\frac{1}{r^2}
       \frac{\partial(r^2\kappa_{rr})}{\partial r} +
     \frac{1}{r}\frac{\partial\kappa_{r\theta}}{\partial\theta}+
     \frac{\cos\theta}{r\sin\theta}\kappa_{r\theta}+
     \frac{1}{r\sin\theta}\frac{\partial\kappa_{r\phi}}{\partial\phi}\right)F_1\right]
   \end{aligned}

The 2nd term one the right is expanded to

.. math::

   \begin{aligned}
     & \frac{\partial^2}{\partial r\partial\theta}\left(\frac{\kappa_{r\theta}}{r}F_1\right) +
     \frac{\partial^2}{\partial\theta^2}\left(\frac{\kappa_{\theta\theta}}{r^2}F_1\right) +
     \frac{\partial^2}{\partial\theta\partial\phi}\left(\frac{\kappa_{\theta\phi}}{r^2\sin\theta}F_1\right) \\
     \nonumber
     & -\frac{\partial}{\partial\theta}\left[\left(\frac{1}{r^2}
       \frac{\partial(r\kappa_{r\theta})}{\partial r}+
     \frac{1}{r^2\sin\theta}\frac{\partial(\sin\theta\kappa_{\theta\theta})}{\partial\theta}+
     \frac{1}{r^2\sin\theta}\frac{\partial\kappa_{\theta\phi}}{\partial\phi}\right)F_1\right]
   \end{aligned}

The 3rd term one the right is expanded to

.. math::

   \begin{aligned}
     & \frac{\partial^2}{\partial r\partial\phi}\left(\frac{\kappa_{r\phi}}{r\sin\theta}F_1\right) +
     \frac{\partial^2}{\partial\theta\partial\phi}\left(\frac{\kappa_{\theta\phi}}{r^2\sin\theta}F_1\right) +
     \frac{\partial^2}{\partial\phi^2}\left(\frac{\kappa_{\phi\phi}}{r^2\sin^2\theta}F_1\right) \\
     \nonumber
     & -\frac{\partial}{\partial\phi}\left[\left(\frac{1}{r^2\sin\theta}
       \frac{\partial(r\kappa_{r\phi})}{\partial r}+
     \frac{1}{r^2\sin\theta}\frac{\partial(\kappa_{\theta\phi})}{\partial\theta}+
     \frac{1}{r^2\sin^2\theta}\frac{\partial\kappa_{\phi\phi}}{\partial\phi}\right)F_1\right]
   \end{aligned}

The final transferred version of Parker transport equation is

.. math::

   \begin{aligned}
     \frac{\partial F_1}{\partial t} = &
     -\frac{\partial}{\partial r}\left[\left(v_r+v_{dr}
       +\frac{1}{r^2}\frac{\partial(r^2\kappa_{rr})}{\partial r} +
     \frac{1}{r}\frac{\partial\kappa_{r\theta}}{\partial\theta}+
     \frac{\cos\theta}{r\sin\theta}\kappa_{r\theta}+
     \frac{1}{r\sin\theta}\frac{\partial\kappa_{r\phi}}{\partial\phi}
     \right)F_1\right] \\\nonumber
     & -\frac{\partial}{\partial\theta}\left[\left(\frac{v_\theta+v_{d\theta}}{r}
       +\frac{1}{r^2}\frac{\partial(r\kappa_{r\theta})}{\partial r}+
     \frac{1}{r^2\sin\theta}\frac{\partial(\sin\theta\kappa_{\theta\theta})}{\partial\theta}+
     \frac{1}{r^2\sin\theta}\frac{\partial\kappa_{\theta\phi}}{\partial\phi}
     \right)F_1\right]\\\nonumber
     & -\frac{\partial}{\partial\phi}\left[\left(\frac{v_\phi+v_{d\phi}}{r\sin\theta}
       +\frac{1}{r^2\sin\theta}\frac{\partial(r\kappa_{r\phi})}{\partial r}+
     \frac{1}{r^2\sin\theta}\frac{\partial(\kappa_{\theta\phi})}{\partial\theta}+
     \frac{1}{r^2\sin^2\theta}\frac{\partial\kappa_{\phi\phi}}{\partial\phi}
     \right)F_1\right]\\\nonumber
     & +\frac{\partial}{\partial p}\left(\frac{p}{3}
     \left(\frac{1}{r^2}\frac{\partial(r^2v_r)}{\partial r} +
     \frac{1}{r\sin\theta}\frac{\partial(\sin\theta v_\theta)}{\partial\theta} +
     \frac{1}{r\sin\theta}\frac{\partial v_\phi}{\partial\phi}\right)F_1\right)\\\nonumber
     & +\frac{\partial^2}{\partial r^2}(\kappa_{rr}F_1) +
     \frac{\partial^2}{\partial r\partial\theta}\left(\frac{\kappa_{r\theta}}{r}F_1\right) +
     \frac{\partial^2}{\partial r\partial\phi}\left(\frac{\kappa_{r\phi}}{r\sin\theta}F_1\right)\\\nonumber
     & +\frac{\partial^2}{\partial r\partial\theta}\left(\frac{\kappa_{r\theta}}{r}F_1\right) +
     \frac{\partial^2}{\partial\theta^2}\left(\frac{\kappa_{\theta\theta}}{r^2}F_1\right) +
     \frac{\partial^2}{\partial\theta\partial\phi}\left(\frac{\kappa_{\theta\phi}}{r^2\sin\theta}F_1\right)\\\nonumber
     & +\frac{\partial^2}{\partial r\partial\phi}\left(\frac{\kappa_{r\phi}}{r\sin\theta}F_1\right) +
     \frac{\partial^2}{\partial\theta\partial\phi}\left(\frac{\kappa_{\theta\phi}}{r^2\sin\theta}F_1\right) +
     \frac{\partial^2}{\partial\phi^2}\left(\frac{\kappa_{\phi\phi}}{r^2\sin^2\theta}F_1\right)\nonumber
   \end{aligned}

This corresponds to a set of SDEs.

.. math::

   \begin{aligned}
     dr & = \left(v_r+v_{dr} +
     \frac{\partial\kappa_{rr}}{\partial r} + \frac{2}{r}\kappa_{rr}+
     \frac{1}{r}\frac{\partial\kappa_{r\theta}}{\partial\theta}+
     \frac{\cos\theta}{r\sin\theta}\kappa_{r\theta}+
     \frac{1}{r\sin\theta}\frac{\partial\kappa_{r\phi}}{\partial\phi}
     \right)dt + [P.dW_t]_r \\
     d\theta & = \left(\frac{v_\theta+v_{d\theta}}{r} +
     \frac{1}{r}\frac{\partial\kappa_{r\theta}}{\partial r} + \frac{\kappa_{r\theta}}{r^2}+
     \frac{1}{r^2}\frac{\partial\kappa_{\theta\theta}}{\partial\theta}+
     \frac{\cos\theta}{r^2\sin\theta}\kappa_{\theta\theta}+
     \frac{1}{r^2\sin\theta}\frac{\partial\kappa_{\theta\phi}}{\partial\phi}
     \right)dt + [P.dW_t]_\theta \\
     d\phi & = \left(\frac{v_\phi+v_{d\phi}}{r\sin\theta}+
     \frac{1}{r\sin\theta}\frac{\partial\kappa_{r\phi}}{\partial r} +
     \frac{\kappa_{r\phi}}{r^2\sin\theta} +
     \frac{1}{r^2\sin\theta}\frac{\partial\kappa_{\theta\phi}}{\partial\theta}+
     \frac{1}{r^2\sin^2\theta}\frac{\partial\kappa_{\phi\phi}}{\partial\phi}
     \right)dt + [P.dW_t]_\phi\\
     dp & = -\frac{p}{3}\left(\frac{\partial v_r}{\partial r}+\frac{2v_r}{r}+
     \frac{1}{r}\frac{\partial v_\theta}{\partial\theta} +
     \frac{\cos\theta}{r\sin\theta}v_\theta +
     \frac{1}{r\sin\theta}\frac{\partial v_\phi}{\partial\phi}\right)
   \end{aligned}

where

.. math::

   \begin{aligned}
     & PP^T =
     \begin{bmatrix}
       2\kappa_{rr} & \dfrac{2\kappa_{r\theta}}{r} & \dfrac{2\kappa_{r\phi}}{r\sin\theta} \\
       \dfrac{2\kappa_{r\theta}}{r} & \dfrac{2\kappa_{\theta\theta}}{r^2} &
       \dfrac{2\kappa_{\theta\phi}}{r^2\sin\theta} \\
       \dfrac{2\kappa_{r\phi}}{r\sin\theta} & \dfrac{2\kappa_{\theta\phi}}{r^2\sin\theta} &
       \dfrac{2\kappa_{\phi\phi}}{r^2\sin^2\theta}
     \end{bmatrix}
   \end{aligned}

According to [Pei10]_, one possibility for
:math:`P` is

.. math::

   \begin{aligned}
     \begin{bmatrix}
       \sqrt{\dfrac{\kappa_{rr}\kappa_{\theta\phi}^2+\kappa_{\theta\theta}\kappa_{r\phi}^2
         +\kappa_{\phi\phi}\kappa_{r\theta}^2-2\kappa_{r\phi}\kappa_{r\theta}\kappa_{\theta\phi}
       -\kappa_{rr}\kappa_{\theta\theta}\kappa_{\phi\phi}}
       {0.5(\kappa_{\theta\phi}^2 - \kappa_{\theta\theta}\kappa_{\phi\phi})}}
       & \dfrac{\kappa_{r\phi}\kappa_{\theta\phi}-\kappa_{r\theta}\kappa_{\phi\phi}}
       {\kappa_{\theta\phi}^2 - \kappa_{\theta\theta}\kappa_{\phi\phi}}
       \sqrt{2\kappa_{\theta\theta}-\dfrac{2\kappa_{\theta\phi}^2}{\kappa_{\phi\phi}}}
       & \dfrac{\sqrt{2}\kappa_{r\phi}}{\sqrt{\kappa_{\phi\phi}}} \\
       0 &
       \dfrac{\sqrt{2\left(\kappa_{\theta\theta}-\kappa_{\theta\phi}^2/\kappa_{\phi\phi}\right)}}{r}
       &
       \dfrac{\kappa_{\theta\phi}}{r}\sqrt{\dfrac{2}{\kappa_{\phi\phi}}} \\
       0 & 0 &
       \dfrac{\sqrt{2\kappa_{\phi\phi}}}{r\sin\theta}
     \end{bmatrix}
   \end{aligned}

For 1D probelms, :math:`F_1=fp^2r^2`, and the corresponding SDEs are

.. math::

   \begin{aligned}
     dr & = \left(v_r + \frac{\partial\kappa_{rr}}{\partial r} +
     \frac{2}{r}\kappa_{rr}\right)dt + \sqrt{2\kappa_{rr}}dW_t \\
     dp & = -\frac{p}{3}\left(\frac{\partial v_r}{\partial r}+\frac{2v_r}{r}\right)
   \end{aligned}

For 2D problems, :math:`F_1=fp^2r^2\sin\theta`, and the corresponding
SDEs are

.. math::

   \begin{aligned}
     dr & = \left(v_r +
     \frac{\partial\kappa_{rr}}{\partial r} + \frac{2}{r}\kappa_{rr}+
     \frac{1}{r}\frac{\partial\kappa_{r\theta}}{\partial\theta}+
     \frac{\cos\theta}{r\sin\theta}\kappa_{r\theta}
     \right)dt + [P.dW_t]_r \\
     d\theta & = \left(\frac{v_\theta}{r} +
     \frac{1}{r}\frac{\partial\kappa_{r\theta}}{\partial r} + \frac{\kappa_{r\theta}}{r^2}+
     \frac{1}{r^2}\frac{\partial\kappa_{\theta\theta}}{\partial\theta}+
     \frac{\cos\theta}{r^2\sin\theta}\kappa_{\theta\theta}
     \right)dt + [P.dW_t]_\theta \\
     dp & = -\frac{p}{3}\left(\frac{\partial v_r}{\partial r}+\frac{2v_r}{r}+
     \frac{1}{r}\frac{\partial v_\theta}{\partial\theta} +
     \frac{\cos\theta}{r\sin\theta}v_\theta\right)
   \end{aligned}

where

.. math::

   \begin{aligned}
     & PP^T =
     \begin{bmatrix}
       2\kappa_{rr} & \dfrac{2\kappa_{r\theta}}{r} \\
       \dfrac{2\kappa_{r\theta}}{r} & \dfrac{2\kappa_{\theta\theta}}{r^2}
     \end{bmatrix}
   \end{aligned}

One possibility for :math:`P` is

.. math::

   \begin{aligned}
     & \begin{bmatrix}
       -\dfrac{Q_{--}\sqrt{-Q_{-+}}}{\sqrt{Q_{--}^2+4b^2}} &
       \dfrac{Q_{+-}\sqrt{Q_{++}}}{\sqrt{Q_{+-}^2+4b^2}} \\
       \dfrac{2b\sqrt{-Q_{-+}}}{\sqrt{Q_{--}^2+4b^2}} &
       \dfrac{2b\sqrt{Q_{++}}}{\sqrt{Q_{+-}^2+4b^2}}
     \end{bmatrix}
   \end{aligned}

where

.. math::

   \begin{aligned}
     Q_{++} &=\sqrt{(a-c)^2+4b^2} + (a + c) \\
     Q_{-+} &=\sqrt{(a-c)^2+4b^2} - (a + c) \\
     Q_{+-} &=\sqrt{(a-c)^2+4b^2} + (a - c) \\
     Q_{--} &=\sqrt{(a-c)^2+4b^2} - (a - c)
   \end{aligned}

where :math:`a=\kappa_{rr}`, :math:`b=\kappa_{r\theta}/r`,
:math:`c=\kappa_{\theta\theta}/r^2`.

.. [Jokipii77] Jokipii, J.R. and Levy, E.H., 1977. Effects of particle drifts on the solar modulation of galactic cosmic rays. The Astrophysical Journal, 213, pp.L85-L88.
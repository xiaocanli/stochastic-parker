3D Model
--------

The relationship
:math:`\sum_\sigma\alpha_\sigma^\mu\alpha_\sigma^\nu = 2\kappa^{\mu\nu}`
is actually a matrix decomposition. We need to decompose
:math:`2\kappa=PP^T`, where
:math:`P=(\boldsymbol{\alpha}_1, \boldsymbol{\alpha}_2, \boldsymbol{\alpha}_3)`.
In a 2D problem, the third component of :math:`\boldsymbol{\alpha}_i` is
essentially 0. In a 3D problem, we need to find all three components of
:math:`\boldsymbol{\alpha}_i`. We need some linear algebra for that.
Every real symmetric matrix can be written in the form
(https://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix#Real_symmetric_matrices)

.. math:: A=Q\Lambda Q^T

where :math:`Q` is an orthogonal matrix whose columns are the
eigenvectors of :math:`A`, and :math:`\Lambda` is a diagonal matrix
whose entries are the eigenvalues of :math:`A`. If the eigenvalues are
non-negative, then the real matrix :math:`P=Q\Lambda^{1/2}`, and

.. math:: A=Q\Lambda^{1/2}\Lambda^{1/2}Q^T = \frac{PP^T}{2}

According to WolframAlpha, the eigenvalue of :math:`\kappa` is
:math:`k_\parallel`, :math:`k_\perp`, and :math:`k_\perp`, and the
corresponding eigenvectors are

.. math::

   \begin{aligned}
     v_1 = \left(\frac{b_x}{b_z}, \frac{b_y}{b_z}, 1\right), \\
     v_2 = \left(-\frac{b_z}{b_x}, 0, 1\right), \\
     v_3 = \left(-\frac{b_y}{b_x}, 1, 0\right).
   \end{aligned}

where :math:`b_x=B_x/B`, :math:`b_y=B_y/B`, and :math:`b_z=B_z/B`.
:math:`v_1`, :math:`v_2`, and :math:`v_3` are not unit vectors, and
:math:`v_2` and :math:`v_3` are not orthogonal to :math:`v_1`, so we
need to re-organize :math:`v_2` and :math:`v_3` and normalize the
vectors.

.. math::

   \begin{aligned}
     v_1 & = \left(b_x, b_y, b_z\right), \\
     v_2 & = \left(-\frac{b_xb_z}{\sqrt{b_x^2+b_y^2}},
     -\frac{b_yb_z}{\sqrt{b_x^2+b_y^2}}, \sqrt{b_x^2+b_y^2}\right),\\
       v_3 & = \left(-\frac{b_y}{\sqrt{b_x^2+b_y^2}}, \frac{b_x}{\sqrt{b_x^2+b_y^2}}, 0\right)
   \end{aligned}

where :math:`v_2` is calculated from the perpendicular component of the
original :math:`v_2` w.r.t. :math:`v_3`. Then,

.. math::

   Q =
     \begin{pmatrix}
       b_x & -b_xb_z/\sqrt{b_x^2+b_y^2} & -b_y/\sqrt{b_x^2+b_y^2}\\
       b_y & -b_yb_z/\sqrt{b_x^2+b_y^2} & b_x/\sqrt{b_x^2+b_y^2}\\
       b_z & \sqrt{b_x^2+b_y^2}         & 0
     \end{pmatrix}

.. math::

   \Lambda =
     \begin{pmatrix}
       \kappa_\parallel & 0 & 0\\
       0 & \kappa_\perp & 0 \\
       0 & 0 & \kappa_\perp
     \end{pmatrix}

.. math::

   P = \sqrt{2}Q\Lambda^{1/2} =
     \begin{pmatrix}
       b_x\sqrt{2\kappa_\parallel} & -b_xb_z\sqrt{2\kappa_\perp}/\sqrt{b_x^2+b_y^2} &
       -b_y\sqrt{2\kappa_\perp}/\sqrt{b_x^2+b_y^2}\\
       b_y\sqrt{2\kappa_\parallel} & -b_yb_z\sqrt{2\kappa_\perp}/\sqrt{b_x^2+b_y^2} &
       b_x\sqrt{2\kappa_\perp}/\sqrt{b_x^2+b_y^2}\\
       b_z\sqrt{2\kappa_\parallel} & \sqrt{b_x^2+b_y^2}\sqrt{2\kappa_\perp}         & 0
     \end{pmatrix}

We can verify that :math:`PP^T=2\kappa`. For 3D simulation, we need to
calculate more terms of the gradient of the diffusion tensor. The
parameters used at particle locations are calculated from :math:`V_x`,
:math:`V_y`, :math:`V_z`, :math:`b_x`, :math:`b_y`, :math:`b_z`,
:math:`\nabla\cdot\boldsymbol{V}`, :math:`\partial_x b_x`,
:math:`\partial_y b_x`, :math:`\partial_z b_x`, :math:`\partial_x b_y`,
:math:`\partial_y b_y`, :math:`\partial_z b_y`, :math:`\partial_x b_z`,
:math:`\partial_y b_z`, :math:`\partial_z b_z`.

.. math::

   \begin{aligned}
     \partial_x\kappa_{xx} & = \partial_x\kappa_\perp +
     \partial_x(\kappa_\parallel-\kappa_\perp)b_x^2 +
     2(\kappa_\parallel-\kappa_\perp)b_x\partial_xb_x, \\
     \partial_y\kappa_{yy} & = \partial_y\kappa_\perp +
     \partial_y(\kappa_\parallel-\kappa_\perp)b_y^2 +
     2(\kappa_\parallel-\kappa_\perp)b_y\partial_yb_y, \\
     \partial_z\kappa_{zz} & = \partial_z\kappa_\perp +
     \partial_z(\kappa_\parallel-\kappa_\perp)b_z^2 +
     2(\kappa_\parallel-\kappa_\perp)b_z\partial_zb_z, \\
     \partial_x\kappa_{xy} & =
     \partial_x(\kappa_\parallel-\kappa_\perp)b_xb_y +
     (\kappa_\parallel-\kappa_\perp)(\partial_xb_xb_y + b_x\partial_xb_y), \\
     \partial_y\kappa_{xy} & =
     \partial_y(\kappa_\parallel-\kappa_\perp)b_xb_y +
     (\kappa_\parallel-\kappa_\perp)(\partial_yb_xb_y + b_x\partial_yb_y), \\
     \partial_x\kappa_{xz} & =
     \partial_x(\kappa_\parallel-\kappa_\perp)b_xb_z +
     (\kappa_\parallel-\kappa_\perp)(\partial_xb_xb_z + b_x\partial_xb_z), \\
     \partial_z\kappa_{xz} & =
     \partial_z(\kappa_\parallel-\kappa_\perp)b_xb_z +
     (\kappa_\parallel-\kappa_\perp)(\partial_zb_xb_z + b_x\partial_zb_z), \\
     \partial_y\kappa_{yz} & =
     \partial_y(\kappa_\parallel-\kappa_\perp)b_yb_z +
     (\kappa_\parallel-\kappa_\perp)(\partial_yb_yb_z + b_y\partial_yb_z), \\
     \partial_z\kappa_{yz} & =
     \partial_z(\kappa_\parallel-\kappa_\perp)b_yb_z +
     (\kappa_\parallel-\kappa_\perp)(\partial_zb_yb_z + b_y\partial_zb_z)
   \end{aligned}

Or we may prefer to use current code structure that calculates
:math:`\partial_x B_x`, :math:`\partial_y B_x`, :math:`\partial_z B_x`,
:math:`\partial_x B_y`, :math:`\partial_y B_y`, :math:`\partial_z B_y`,
:math:`\partial_x B_z`, :math:`\partial_y B_z`, :math:`\partial_z B_z`.
Then, the derivatives are calculated as

.. math::

   \begin{aligned}
     \partial_xB & = b_x\partial_xB_x + b_y\partial_xB_y + b_z\partial_xB_z, \\
     \partial_yB & = b_x\partial_yB_x + b_y\partial_yB_y + b_z\partial_yB_z, \\
     \partial_zB & = b_x\partial_zB_x + b_y\partial_zB_y + b_z\partial_zB_z, \\
     \partial_x\kappa_{xx} & = \partial_x\kappa_\perp +
     \partial_x(\kappa_\parallel-\kappa_\perp)b_x^2 +
     2(\kappa_\parallel-\kappa_\perp)(B_xB\partial_xB_x - B_x^2\partial_x B)/B^3, \\
     \partial_y\kappa_{yy} & = \partial_y\kappa_\perp +
     \partial_y(\kappa_\parallel-\kappa_\perp)b_y^2 +
     2(\kappa_\parallel-\kappa_\perp)(B_yB\partial_yB_y - B_y^2\partial_y B)/B^3, \\
     \partial_z\kappa_{zz} & = \partial_z\kappa_\perp +
     \partial_z(\kappa_\parallel-\kappa_\perp)b_z^2 +
     2(\kappa_\parallel-\kappa_\perp)(B_zB\partial_zB_z - B_z^2\partial_z B)/B^3, \\
     \partial_x\kappa_{xy} & = \partial_x(\kappa_\parallel-\kappa_\perp)b_xb_y +
     (\kappa_\parallel-\kappa_\perp)[(B_y\partial_xB_x + B_x\partial_xB_y)B -
     2B_xB_y\partial_xB] / B^3, \\
     \partial_y\kappa_{xy} & = \partial_y(\kappa_\parallel-\kappa_\perp)b_xb_y +
     (\kappa_\parallel-\kappa_\perp)[(B_y\partial_yB_x + B_x\partial_yB_y)B -
     2B_xB_y\partial_yB] / B^3, \\
     \partial_x\kappa_{xz} & = \partial_x(\kappa_\parallel-\kappa_\perp)b_xb_z +
     (\kappa_\parallel-\kappa_\perp)[(B_z\partial_xB_x + B_x\partial_xB_z)B -
     2B_xB_z\partial_xB] / B^3, \\
     \partial_z\kappa_{xz} & = \partial_z(\kappa_\parallel-\kappa_\perp)b_xb_z +
     (\kappa_\parallel-\kappa_\perp)[(B_z\partial_zB_x + B_x\partial_zB_z)B -
     2B_xB_z\partial_zB] / B^3, \\
     \partial_y\kappa_{yz} & = \partial_y(\kappa_\parallel-\kappa_\perp)b_yb_z +
     (\kappa_\parallel-\kappa_\perp)[(B_z\partial_yB_y + B_y\partial_yB_z)B -
     2B_yB_z\partial_yB] / B^3, \\
     \partial_z\kappa_{yz} & = \partial_z(\kappa_\parallel-\kappa_\perp)b_yb_z +
     (\kappa_\parallel-\kappa_\perp)[(B_z\partial_zB_y + B_y\partial_zB_z)B -
     2B_yB_z\partial_zB] / B^3.
   \end{aligned}

Particle drift velocity
^^^^^^^^^^^^^^^^^^^^^^^

In the 3D model, we need the drift velocity, which is given by

.. math::

   \begin{aligned}
     & \boldsymbol{V}_d = \frac{pcw}{3q}\nabla\times\left(\frac{\boldsymbol{B}}{B^2}\right)
     = \frac{1}{3q}\frac{p^2c}{\sqrt{p^2+m^2c^2}}
     \left(\frac{1}{B^2}\nabla\times\boldsymbol{B} -
     \frac{2}{B^3}\nabla B\times\boldsymbol{B}\right) \\
     & \nabla\times\boldsymbol{B} =
     (\partial_y B_z - \partial_z B_y)\hat{i} +
     (\partial_z B_x - \partial_x B_z)\hat{j} +
     (\partial_x B_y - \partial_y B_x)\hat{k} \\
     & \nabla B\times\boldsymbol{B} =
     (B_z\partial_yB - B_y\partial_zB)\hat{i} +
     (B_x\partial_zB - B_z\partial_xB)\hat{j} +
     (B_y\partial_xB - B_x\partial_yB)\hat{k}
   \end{aligned}

where :math:`p=\gamma m v` is particle momentum, :math:`c` is the speed
of light, :math:`w=v/c` is the normalized particle speed, and :math:`q`
is particle charge. Using normalized quantities, we have

.. math::

   \begin{aligned}
     \tilde{\boldsymbol{V}}_d & = \frac{1}{v_A}\frac{1}{3\tilde{q}e}\frac{\tilde{p}^2p_0^2c}{\sqrt{\tilde{p}^2p_0^2+m^2c^2}}\frac{1}{B_0L_0}
     \left(\frac{1}{\tilde{B}^2}\tilde{\nabla}\times\tilde{\boldsymbol{B}} -
     \frac{2}{\tilde{B}^3}\tilde{\nabla}\tilde{B}\times\tilde{\boldsymbol{B}}\right) \\
     & = \frac{1}{\sqrt{d_1^2\tilde{p}^{-2}+d_2^2\tilde{p}^{-4}}}
     \frac{1}{3\tilde{q}}\left(\frac{1}{\tilde{B}^2}\tilde{\nabla}\times\tilde{\boldsymbol{B}} -
     \frac{2}{\tilde{B}^3}\tilde{\nabla}\tilde{B}\times\tilde{\boldsymbol{B}}\right)
   \end{aligned}

where :math:`\tilde{\boldsymbol{V}}_d=\boldsymbol{V}_d/v_A`,
:math:`\tilde{q}=q/e`, :math:`\tilde{\nabla}=L_0\nabla`,
:math:`\tilde{\boldsymbol{B}}=\boldsymbol{B}/B_0`,
:math:`\tilde{p}=p/p_0`, :math:`d_1=eB_0v_AL_0/(p_0c)`, and
:math:`d_2=emB_0v_AL_0/p_0^2`. Note that in the code, :math:`\tilde{p}`
will be re-normalized. For example, :math:`\tilde{p}_0=1` might
correspond to :math:`\tilde{p}_{n0}=0.1` in simulations. The
re-normalized numerical momentum
:math:`\tilde{p}_n=\tilde{p}\tilde{p}_{n0}`. Thus,
:math:`\tilde{p} = \tilde{p}_n/\tilde{p}_{n0}` in simulations, and we
need provide :math:`d_1` and :math:`d_2` based on the normalization.

.. note::
  The velocity normalization :math:`v_A` should be changed to :math:`v_0` if :math:`v_0\neq v_A`
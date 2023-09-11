Algorithm & Theory
==================

.. autosummary::
   :toctree: generated

.. role:: raw-latex(raw)
   :format: latex
..

1D & 2D model
-------------

Parker’s transport equation

.. math::
   :name: equ_parker

   \frac{\partial f}{\partial t} + (\boldsymbol{v}+\boldsymbol{v}_d)\cdot\nabla f
     - \frac{1}{3}\nabla\cdot\boldsymbol{v}\frac{\partial f}{\partial\ln p}
     = \nabla\cdot(\boldsymbol{\kappa}\nabla f) + Q,

where :math:`f(x_i, p, t)` is the particle distribution function as a
function of the particle position :math:`x_i`, momentum :math:`p`
(isotropic momentum assumed), and time :math:`t`;
:math:`\boldsymbol{\kappa}` is the spatial diffusion coefficient tensor,
:math:`\boldsymbol{v}` is the bulk plasma velocity, and :math:`Q` is the
source. Note that the particle drift :math:`\boldsymbol{v}_d` is out of
the simulation plane and is not considered here. The diffusion
coefficient tensor is given by

.. math::

   \kappa_{ij} = \kappa_\perp\delta_{ij} -
     \frac{(\kappa_\perp-\kappa_\parallel)B_iB_j}{B^2},

where :math:`\kappa_\parallel` and :math:`\kappa_\perp` are the parallel
and perpendicular diffusion coefficients. Here :math:`\kappa_\parallel`
can be calculated from the quasi-linear
theory [8]_. Assuming that
magnetic turbulence is well developed and has an isotropic power
spectrum :math:`P\sim k^{-5/3}`, the resulting
:math:`\kappa_\parallel\sim p^{4/3}` when the particle gyroradius is
much smaller than the correlation length of turbulence. In particular,
we use the following expression for
:math:`\kappa_\parallel` [7]_,

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
variance of turbulence. Reference [7]_ gave
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
fluctuations [1]_) affect particle
parallel transport [6]_. The
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
energy [7]_. There is also
observational evidence suggesting that
:math:`\kappa_\perp/\kappa_\parallel` can be quite
large [3]_ [21]_.

The Parker transport equation can be solved by integrating the
stochastic differential equation corresponding to the Fokker–Planck form
of the transport
equation [20]_ [5]_ [12]_ [10]_.
Neglecting the source term :math:`Q` in
equ_parker_ and assuming :math:`F=fp^2`,

.. math::

   \begin{aligned}
     \frac{\partial F}{\partial t}
     & = -\nabla\cdot\left[(\nabla\cdot\boldsymbol{\kappa}+\boldsymbol{v})F\right] +
     \frac{\partial}{\partial p} \left[\frac{p}{3}\nabla\cdot\boldsymbol{v} F\right] +
     \nabla\cdot(\nabla\cdot(\boldsymbol{\kappa}F)),
   \end{aligned}

which is equivalent to a system of stochastic differential equations
(SDEs) of the Ito type. In a 1D model,

.. math::

   dX = \left(\frac{\partial\kappa}{\partial x} + v_x\right)ds +
     \sqrt{2\kappa} dW_\sigma(s),\quad
     dp=-\frac{p}{3}\frac{\partial v_x}{\partial x}ds,

where :math:`dW` is the normalized distributed random number with mean
zero and variance :math:`\sqrt{\Delta t}`, and :math:`\Delta t` is the
time step for stochastic integration. This corresponds to a Wiener
process. Numerical approximation is often used for the Wiener process to
replace the normal distribution. We use a uniform distribution in
:math:`[-\sqrt{3}, \sqrt{3}]` in the code.

In a 2D model,

.. math::

   dX = (\nabla\cdot\boldsymbol{\kappa} + \boldsymbol{v})ds +
     \sum_\sigma\boldsymbol{\alpha}_\sigma dW_\sigma(s),\quad
     dp=-\frac{p}{3}(\nabla\cdot\boldsymbol{v})ds,

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
     \end{pmatrix}, \quad

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
The detailed expressions depend on the diffusion model to choose. Using
equ_parker_,

.. math::

   \begin{aligned}
     \frac{\partial\kappa}{\partial x}\sim\kappa\left(
     \frac{2}{3L_c}\frac{\partial L_c}{\partial x} -
     \frac{1}{3B}\frac{\partial B}{\partial x} -
     \frac{1}{\sigma^2}\frac{\partial(\sigma^2)}{\partial x}
     \right)
   \end{aligned}

In this work, we use a derivative-free Milstein
method [2]_ to solve the stochastic
differential equation. It is different from the usual method due to one
more term, which makes it become a higher-order method:

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
computation. For a 1D problem, the particle moves a distance satisfying
:math:`l_x^2=\text{max}\left(\left<\Delta x\right>^2,
\left<\Delta x^2\right>\right)` [16]_,
where

.. math::

   \begin{aligned}
     \left<\Delta x\right> = \left(v_x + \frac{d\kappa(x)}{dx}\right)\Delta t,
     \quad \left<\Delta x^2\right> = 2\kappa(x)\Delta t,
   \end{aligned}

and :math:`l_x` should be much smaller than the spatial variation scale
of the fields. In this work, we assume
:math:`\left<\Delta x\right>^2 < \left<\Delta x^2\right>` and choose
:math:`\Delta t` so that :math:`l_x\ll\delta_x`, where :math:`\delta_x`
is the grid size. For our 2D problems, we choose the following criteria
to determine the time step:

.. math::

   \begin{aligned}
     \Delta t_x & = \text{min}\left[\frac{\delta x}{80|v_x + \partial_x\kappa_{xx} +
     \partial_y\kappa_{xy}|},
     \frac{\left(\sqrt{2\kappa_\perp} + \sqrt{2(\kappa_\parallel - \kappa_\perp)}|B_x/B|\right)^2}
     {|v_x + \partial_x\kappa_{xx} + \partial_y\kappa_{xy}|^2}\right], \\
     \Delta t_y & = \text{min}\left[\frac{\delta y}{80|v_y + \partial_y\kappa_{yy} +
     \partial_x\kappa_{xy}|},
     \frac{\left(\sqrt{2\kappa_\perp} + \sqrt{2(\kappa_\parallel - \kappa_\perp)}|B_y/B|\right)^2}
     {|v_y + \partial_y\kappa_{yy} + \partial_x\kappa_{xy}|^2}\right],\\
     \Delta t & = \text{min}(\Delta t_x, \Delta t_y).
   \end{aligned}

3D model
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
parameters used at particle locations are calculated from :math:`v_x`,
:math:`v_y`, :math:`v_z`, :math:`b_x`, :math:`b_y`, :math:`b_z`,
:math:`\nabla\cdot\boldsymbol{v}`, :math:`\partial_x b_x`,
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

In 3D, we need the drift velocity, which is given by

.. math::

   \begin{aligned}
     & \boldsymbol{v}_d = \frac{pcw}{3q}\nabla\times\left(\frac{\boldsymbol{B}}{B^2}\right)
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
     \tilde{\boldsymbol{v}}_d & = \frac{1}{v_A}\frac{1}{3\tilde{q}e}\frac{\tilde{p}^2p_0^2c}{\sqrt{\tilde{p}^2p_0^2+m^2c^2}}\frac{1}{B_0L_0}
     \left(\frac{1}{\tilde{B}^2}\tilde{\nabla}\times\tilde{\boldsymbol{B}} -
     \frac{2}{\tilde{B}^3}\tilde{\nabla}\tilde{B}\times\tilde{\boldsymbol{B}}\right) \\
     & = \frac{1}{\sqrt{d_1^2\tilde{p}^{-2}+d_2^2\tilde{p}^{-4}}}
     \frac{1}{3\tilde{q}}\left(\frac{1}{\tilde{B}^2}\tilde{\nabla}\times\tilde{\boldsymbol{B}} -
     \frac{2}{\tilde{B}^3}\tilde{\nabla}\tilde{B}\times\tilde{\boldsymbol{B}}\right)
   \end{aligned}

where :math:`\tilde{\boldsymbol{v}}_d=\boldsymbol{v}_d/v_A`,
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

Momentum Diffusion
------------------

We can include an momentum diffusion term to the right side of
equ_parker_.

.. math::
   :name: equ_parker_2nd

   \frac{\partial f}{\partial t} + (\boldsymbol{v}+\boldsymbol{v}_d)\cdot\nabla f
     - \frac{1}{3}\nabla\cdot\boldsymbol{v}\frac{\partial f}{\partial\ln p}
     = \nabla\cdot(\boldsymbol{\kappa}\nabla f) +
     \frac{1}{p^2}\frac{\partial}{\partial p}
     \left(p^2D_{pp}\frac{\partial f}{\partial p}\right) + Q,

Neglecting the source term :math:`Q` in
equ_parker_2nd_ and assuming :math:`F=fp^2`,

.. math::

   \begin{aligned}
     \frac{\partial F}{\partial t} =
     & -\nabla\cdot\left[(\nabla\cdot\boldsymbol{\kappa}+\boldsymbol{v}+\boldsymbol{v}_d)F\right] +
     \nabla\cdot(\nabla\cdot(\boldsymbol{\kappa}F)) + \nonumber \\
     & \frac{\partial}{\partial p} \left[\left(\frac{p}{3}\nabla\cdot\boldsymbol{v} -
     \frac{\partial D_{pp}}{\partial p} - \frac{2D_{pp}}{p}\right) F\right] +
     \frac{\partial(D_{pp}F)}{\partial p^2}.
   \end{aligned}

which is equivalent to a system of stochastic differential equations
(SDEs) of the Ito type,

.. math::

   \begin{aligned}
     dX & = (\nabla\cdot\boldsymbol{\kappa} + \boldsymbol{v} + \boldsymbol{v}_d)ds +
     \sum_\sigma\boldsymbol{\alpha}_\sigma dW_\sigma(s) \\
     dp & = \left(-\frac{p}{3}\nabla\cdot\boldsymbol{v} +
     \frac{\partial D_{pp}}{\partial p} + \frac{2D_{pp}}{p}\right)ds +
     \sqrt{2D_{pp}}dW(s)
   \end{aligned}

where
:math:`\sum_\sigma\alpha_\sigma^\mu\alpha_\sigma^\nu = 2\kappa^{\mu\nu}`,
:math:`dW` is the normalized distributed random number with mean zero
and variance :math:`\sqrt{\Delta t}`, and :math:`\Delta t` is the time
step for stochastic integration. This corresponds to a Wiener process.
Numerical approximation is often used for the Wiener process to replace
the normal distribution. We use a uniform distribution in
:math:`[-\sqrt{3}, \sqrt{3}]` in the code. For a 2D problem,
reference [15]_ shows that for forward and
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
models [13]_ [14]_ [11]_. The corresponding SDE is

.. math::

   \begin{aligned}
     dp & = \left(-\frac{p}{3}\nabla\cdot\boldsymbol{v} + \frac{4pv_A^2}{9\kappa_\parallel}\right)ds +
     \sqrt{\frac{2p^2v_A^2}{9\kappa_\parallel}}dW(s), \text{if $\kappa_\parallel$ is independent of $p$}\\
     dp & = \left(-\frac{p}{3}\nabla\cdot\boldsymbol{v} + \frac{8pv_A^2}{27\kappa_\parallel}\right)ds +
     \sqrt{\frac{2p^2v_A^2}{9\kappa_\parallel}}dW(s), \text{if $\kappa_\parallel\sim p^{4/3}$}
   \end{aligned}

which are normalized to

.. math::

   \begin{aligned}
     d\tilde{p}_n & = \left(-\frac{\tilde{p}_n}{3}\tilde{\nabla}\cdot\tilde{\boldsymbol{v}} + \frac{4\tilde{p}_n\tilde{v}_A^2}{9\tilde{\kappa}_\parallel}\right)d\tilde{s} + \tilde{p}_n\tilde{v}_A\sqrt{\frac{2}{9\tilde{\kappa}_\parallel}}dW(\tilde{s}), \text{if $\kappa_\parallel$ is independent of $p$}\\
     d\tilde{p}_n & = \left(-\frac{\tilde{p}_n}{3}\tilde{\nabla}\cdot\tilde{\boldsymbol{v}} + \frac{8\tilde{p}_n\tilde{v}_A^2}{27\tilde{\kappa}_\parallel}\right)d\tilde{s} + \tilde{p}_n\tilde{v}_A\sqrt{\frac{2}{9\tilde{\kappa}_\parallel}}dW(\tilde{s}), \text{if $\kappa_\parallel\sim p^{4/3}$}
   \end{aligned}

where :math:`\tilde{p}_n=\tilde{p}\tilde{p}_{n0}=p\tilde{p}_{n0}/p_0`,
where is :math:`\tilde{p}_{n0}` is the numerical value for particles
with :math:`p_0` in the code (e.g., 0.1 as often used),
:math:`\tilde{\nabla}=L_0\nabla`,
:math:`\tilde{\boldsymbol{v}}=\boldsymbol{v}/v_{A0}`,
:math:`\tilde{v}_A=\tilde{v}_{A0}`,
:math:`\tilde{\kappa}_\parallel=\kappa_\parallel/\kappa_0`,
:math:`\kappa_0=L_0v_{A0}`, :math:`\tilde{s}=s/t_0`, and
:math:`t_0=L_0/v_{A0}`. These are all given in the code.

For isotropic particle distributions, the flow shear introduces another
momentum diffusion term. If there is no average magnetic
field [4]_.

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
scattering. According to [17]_, :math:`\tau`
is related particle diffusion coefficient
:math:`\kappa_\parallel=v^2\tau/3`. The corresponding SDE is

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
:math:`\tau\sim\tau_0(p_0/p)^2` [4]_,

.. math::

   \begin{aligned}
     d\tilde{p}_n = \left(-\frac{\tilde{p}_n}{3}\tilde{\nabla}\cdot\tilde{\boldsymbol{v}} + \frac{2\tilde{\Gamma}\tilde{\tau}_0\tilde{p}_{n0}^2}{\tilde{p}_n}\right)d\tilde{s} + \sqrt{2\tilde{\Gamma}\tilde{\tau}_0\tilde{p}_{n0}^2}dW(\tilde{s})
   \end{aligned}

For
:math:`\tau\sim\tau_0(p_0/p)^{2/3}` [7]_,

.. math::

   \begin{aligned}
     d\tilde{p}_n & = \left(-\frac{\tilde{p}_n}{3}\tilde{\nabla}\cdot\tilde{\boldsymbol{v}} + \frac{10}{3}\tilde{\Gamma}\tilde{\tau}_0\tilde{p}_{n}^{1/3}\tilde{p}_{n0}^{2/3}\right)d\tilde{s} + \sqrt{2\tilde{\Gamma}\tilde{\tau}_0\tilde{p}_n^{4/3}\tilde{p}_{n0}^{2/3}}dW(\tilde{s}) \\
     \tau_0 & = 3\kappa_{\parallel 0} / v_0^2
   \end{aligned}

If there is an average magnetic field, the equation is more complicated
(see [18]_ [19]_).

Spherical Coordinates
---------------------

In spherical coordinates, the drift velocity

.. math::

   \begin{aligned}
     & \boldsymbol{v}_d = \frac{pcw}{3q}\nabla\times\left(\frac{\boldsymbol{B}}{B^2}\right)
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
     & (\nabla B\times\vect{B})_r = (\nabla B)_\theta B_\phi - (\nabla B)_\phi B_\theta \\
     & (\nabla B\times\vect{B})_\theta = (\nabla B)_\phi B_r - (\nabla B)_r B_\phi \\
     & (\nabla B\times\vect{B})_\phi = (\nabla B)_r B_\theta - (\nabla B)_\theta B_r
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
coordinates. Since we don’t have cross diffusion terms (spatial and
momentum), we can ignore the momentum diffusion for now.

.. math::

   \begin{aligned}
     \frac{\partial F}{\partial t}
     & = -(\boldsymbol{v}+\boldsymbol{v}_d)\cdot\nabla F
     - (\nabla\cdot\boldsymbol{v})F
     + \frac{\partial}{\partial p}\left[\frac{p}{3}(\nabla\cdot\boldsymbol{v})F\right]
     + \nabla\cdot(\boldsymbol{\kappa}\cdot\nabla F)
   \end{aligned}

where :math:`F=fp^2`. Since :math:`\nabla\cdot\boldsymbol{v}_d=0`, we
can add one more term :math:`-(\nabla\cdot\boldsymbol{v}_d)F` to the
right. Then,

.. math::

   \begin{aligned}
     \frac{\partial F}{\partial t}
     & = -\nabla\cdot((\boldsymbol{v}+\boldsymbol{v}_d)F)
     + \frac{\partial}{\partial p}\left[\frac{p}{3}(\nabla\cdot\boldsymbol{v})F\right]
     + \nabla\cdot(\boldsymbol{\kappa}\cdot\nabla F)
   \end{aligned}

Taking :math:`\boldsymbol{V}=\boldsymbol{v}+\boldsymbol{v}_d`,

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
:math:`F_1=F\sin\theta r^2` [9]_ [12]_.
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

According to [12]_, one possibility for
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

.. [1] Bieber, J.W., Wanner, W. and Matthaeus, W.H., 1996. Dominant two‐dimensional solar wind turbulence with implications for cosmic ray transport. Journal of Geophysical Research: Space Physics, 101(A2), pp.2511-2522.
.. [2] Burrage, K., Burrage, P.M. and Tian, T., 2004. Numerical methods for strong solutions of stochastic differential equations: an overview. Proceedings of the Royal Society of London. Series A: Mathematical, Physical and Engineering Sciences, 460(2041), pp.373-402.
.. [3] Dwyer, J.R., Mason, G.M., Mazur, J.E., Jokipii, J.R., Von Rosenvinge, T.T. and Lepping, R.P., 1997. Perpendicular transport of low-energy corotating interaction region-associated nuclei. The Astrophysical Journal, 490(1), p.L115.
.. [4] Earl, J.A., Jokipii, J.R. and Morfill, G., 1988. Cosmic-ray viscosity. The Astrophysical Journal, 331, pp.L91-L94.
.. [5] Florinski, V. and Pogorelov, N.V., 2009. Four-dimensional transport of galactic cosmic rays in the outer heliosphere and heliosheath. The Astrophysical Journal, 701(1), p.642.
.. [6] Florinski, V., Zank, G.P. and Pogorelov, N.V., 2003. Galactic cosmic ray transport in the global heliosphere. Journal of Geophysical Research: Space Physics, 108(A6).
.. [7] Giacalone, J. and Jokipii, J.R., 1999. The transport of cosmic rays across a turbulent magnetic field. The Astrophysical Journal, 520(1), p.204.
.. [8] Jokipii, J.R., 1971. Propagation of cosmic rays in the solar wind. Reviews of Geophysics, 9(1), pp.27-87.
.. [9] Jokipii, J.R. and Levy, E.H., 1977. Effects of particle drifts on the solar modulation of galactic cosmic rays. The Astrophysical Journal, 213, pp.L85-L88.
.. [10] Kong, X., Guo, F., Giacalone, J., Li, H. and Chen, Y., 2017. The acceleration of high-energy protons at coronal shocks: the effect of large-scale streamer-like magnetic field structures. The Astrophysical Journal, 851(1), p.38.
.. [11] Le Roux, J.A. and Webb, G.M., 2007. Nonlinear cosmic-ray diffusive transport in combined two-dimensional and slab magnetohydrodynamic turbulence: a BGK-Boltzmann approach. The Astrophysical Journal, 667(2), p.930.
.. [12] Pei, C., Bieber, J. W., Burger, R. A., & Clem, J. 2010, Journal of Geophysical Research (Space Physics), 115, A12107
.. [13] Schlickeiser, R., 1989. Cosmic-ray transport and acceleration. I-Derivation of the kinetic equation and application to cosmic rays in static cold media. II-Cosmic rays in moving cold media with application to diffusive shock wave acceleration. The Astrophysical Journal, 336, pp.243-293.
.. [14] Schlickeiser, R. and Miller, J.A., 1998. Quasi-linear theory of cosmic ray transport and acceleration: the role of oblique magnetohydrodynamic waves and transit-time damping. The Astrophysical Journal, 492(1), p.352.
.. [15] Skilling, J., 1975. Cosmic Ray Streaming—II effect of particles on alfvén waves. Monthly Notices of the Royal Astronomical Society, 173(2), pp.245-254.
.. [16] Strauss, R. and Effenberger, F., 2017. A hitch-hiker’s guide to stochastic differential equations. Space Science Reviews, 212(1), pp.151-192.
.. [17] Webb, G. M., Barghouty, A. F., Hu, Q., & le Roux, J. A. 2018, The Astrophysical Journal, 855, 31
.. [18] Williams, L.L. and Jokipii, J.R., 1991. Viscosity and inertia in cosmic-ray transport-Effects of an average magnetic field. The Astrophysical Journal, 371, pp.639-647.
.. [19] Williams, L.L., Schwadron, N., Jokipii, J.R. and Gombosi, T.I., 1993. A unified transport equation for both cosmic rays and thermal particles. The Astrophysical Journal, 405, pp.L79-L81.
.. [20] Zhang, M., 1999. A Markov stochastic process theory of cosmic-ray modulation. The Astrophysical Journal, 513(1), p.409.
.. [21] Zhang, M., Jokipii, J.R. and McKibben, R.B., 2003. Perpendicular transport of solar energetic particles in heliospheric magnetic fields. The Astrophysical Journal, 595(1), p.493.

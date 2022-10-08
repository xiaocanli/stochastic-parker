Overview
========

The code was originally from Fan Guo [1]_. It has been updated over the years. See recent papers using this code at :doc:`publications`. The code solves Parker's transport equation.

.. math::
  \frac{\partial f}{\partial t} + (\boldsymbol{V}+\boldsymbol{V}_d)\cdot\nabla f
  - \frac{1}{3}\nabla\cdot\boldsymbol{V}\frac{\partial f}{\partial\ln p}
  = \nabla\cdot(\boldsymbol{\kappa}\nabla f) + Q,

where :math:`f(x_i, p, t)` is the particle distribution function as a function of the particle position :math:`x_i`,  momentum :math:`p` (isotropic momentum assumed), and time :math:`t`; :math:`\boldsymbol{\kappa}` is the spatial diffusion coefficient tensor, :math:`\boldsymbol{V}` is the bulk plasma velocity, :math:`\boldsymbol{V}_d` is the particle drift, and :math:`Q` is the source. In general, :math:`\boldsymbol{V}` can be obtained directly from the MHD simulations, :math:`\boldsymbol{V}_d` and :math:`\boldsymbol{\kappa}` both depend on the vector magnetic field, and :math:`Q` could depend on many plasma properties (e.g., number density, current density, and compression). :math:`\boldsymbol{\kappa}` also depends on the turbulence properties (e.g., amplitude, spectral slope, and anisotropy). See :doc:`theory` for more details about the algorithms.

.. [1] https://ui.adsabs.harvard.edu/abs/2010ApJ...725..128G/abstract
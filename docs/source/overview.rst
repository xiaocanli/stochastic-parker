Overview
========

The code was originally from Fan Guo [1]_. It has been updated over the years. See recent papers using this code at :doc:`publications`. The code solves Parker's transport equation

.. math::
  \frac{\partial f}{\partial t} + (\boldsymbol{v}+\boldsymbol{v}_d)\cdot\nabla f
  - \frac{1}{3}\nabla\cdot\boldsymbol{v}\frac{\partial f}{\partial\ln p}
  = \nabla\cdot(\boldsymbol{\kappa}\nabla f) + Q,

where :math:`f(x_i, p, t)` is the particle distribution function as a function of the particle position :math:`x_i`,  momentum :math:`p` (isotropic momentum assumed), and time :math:`t`; :math:`\boldsymbol{\kappa}` is the spatial diffusion coefficient tensor, :math:`\boldsymbol{v}` is the bulk plasma velocity, :math:`\boldsymbol{v}_d` is the particle drift, and :math:`Q` is the source. In general, :math:`\boldsymbol{v}` can be obtained directly from the MHD simulations, :math:`\boldsymbol{v}_d` and :math:`\boldsymbol{\kappa}` both depend on the vector magnetic field, and :math:`Q` could depend on many plasma properties (e.g., number density, current density, and compression). :math:`\boldsymbol{\kappa}` also depends on the turbulence properties (e.g., amplitude, spectral slope, and anisotropy). See :doc:`theory` for more details about the transport theories.

The code solves the equation using the stochastic integration method since the Fokker-Planck form of the transport equation is equivalent to a set of the stochastic differential equations (SDEs) of Ito type [2]_. The SDEs are solved using a large number of pseudo particles sampling :math:`F=fp^2`. See :doc:`theory` for more details about the algorithms.

The simulations will produce spatially and temporally dependent energy spectra and maps of energetic particles. See :doc:`usage` for how to get these from MHD simulation outputs.

.. [1] https://ui.adsabs.harvard.edu/abs/2010ApJ...725..128G/abstract
.. [2] https://en.wikipedia.org/wiki/Fokker%E2%80%93Planck_equation
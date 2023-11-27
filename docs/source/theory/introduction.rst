Introduction
============

The GPAT code was initially developed by Fan Guo [1]_. It has been updated over the years. Recent papers using this code can be found at :doc:`../publications`. The code solves particle transport equations, including Parker's transport equation and focused transport equation, using the stochastic integration method. It is supposed to be used for studying particle acceleration and transport at scales much larger than the kinetic scales.

The code solves the equation using the stochastic integration method since the Fokker-Planck form of the transport equation is equivalent to a set of the stochastic differential equations (SDEs) of Ito type [2]_. The SDEs are solved using a large number of pseudo particles.

It uses MHD simulation data (magnetic field, velocity, density, etc.) as background to evolve the transport equations. It has been used in studying particle acceleration and transport in magnetic reconnection, solar flares, and solar eruption regions.

The simulations will produce spatially and temporally dependent energy spectra and maps of energetic particles. See :doc:`../usage` for how to get these from MHD simulation outputs.

.. [1] https://ui.adsabs.harvard.edu/abs/2010ApJ...725..128G/abstract
.. [2] https://en.wikipedia.org/wiki/Fokker%E2%80%93Planck_equation
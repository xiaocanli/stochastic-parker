Welcome to GPAT's documentation!
=============================================

**GPAT** (Global Particle Acceleration and Transport)  is a program for solving particle transport equations. It solves either Parker's transport equation for nearly isotropic particle distributions or the focused transport equations for anisotropic particle distributions. The code uses the stochastic integration method to solve these equations. For details on the algorithms implemented, please check out the theory section.

Contact us
----------

If you're beginning to use GPAT or have any questions, feel free to drop by the `discussions page <https://github.com/xiaocanli/stochastic-parker/discussions>`__. For bug reports or to request new features, you can also open a new `issue <https://github.com/xiaocanli/stochastic-parker/issues>`__.

.. raw:: html

   <style>
   /* front page: hide chapter titles
    * needed for consistent HTML-PDF-EPUB chapters
    */
   section#install,
   section#usage,
   section#theory,
   section#development,
   section#epilogue {
       display:none;
   }
   </style>

Install
------------
.. toctree::
   :caption: INSTALLATION
   :maxdepth: 1
   :hidden:

   install

Usage
------------
.. toctree::
   :caption: USAGE
   :maxdepth: 1
   :hidden:

   usage
   input_parameters

Theory
------------
.. toctree::
   :caption: THEORY
   :maxdepth: 1
   :hidden:

   theory/introduction
   theory/parker_transport
   theory/focused_transport

Development
------------
.. toctree::
   :caption: DEVELOPMENT
   :maxdepth: 1
   :hidden:

   development/implementation

Epilogue
------------
.. toctree::
   :caption: EPILOGUE
   :maxdepth: 1
   :hidden:

   publications
   acknowledgments
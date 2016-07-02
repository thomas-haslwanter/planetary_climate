Principles of Planetary Climate - Documentation
===============================================

*Python 3* utilities for the book `Principles of Planetary Climate <https://geosci.uchicago.edu/~rtp1/PrinciplesPlanetaryClimate/>`_, by Raymond T. Pierrehumbert


It is hosted under https://github.com/thomas-haslwanter/planetary_climate, and contains the following folders:

*ClimateUtilities*  Analysis routines for climate simulations
    - Written for stand-alone analysis

*cu_sp (ClimateUtilities for SciPy)* 
    - Requires the Python packages *scipy* and *pandas*


Installation
------------

Just copy the files to your computer, and choose with which set of tools
(*ClimateUtilities* or *cu_sp*) you want to work with. Then go into that
directory, and type on the commandline

    >>> python setpath.py

Dependencies
^^^^^^^^^^^^
numpy, scipy, matplotlib, pandas


Modules
-------

.. toctree::
   :maxdepth: 2

   gases
   math_demos
   phys
   planets
   satvp
   setpath


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. note::
    | *Author:*     Thomas Haslwanter
    | *Version:*    0.1.0
    | *Date:*       June 2016
    | *email:*      thomas.haslwanter@fh-linz.at
    | *Copyright (c):*      2016, Thomas Haslwanter. All rights reserved.
    | *Licence:*    This work is licensed under the `BSD 2-Clause License <http://opensource.org/licenses/BSD-2-Clause>`_


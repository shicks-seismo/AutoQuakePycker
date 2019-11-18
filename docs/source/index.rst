.. AutoQuakePycker documentation master file, created by
   sphinx-quickstart on Fri Nov 15 14:49:15 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

AutoQuakePycker
===========================================

AutoQuakePycker is a Python package for the automatic and iterative picking and relocation
of local seismic events. 
AutoQuakePycker takes a first-guess catalogue of earthquake locations (I recommend the very
easy-to-use and comprehnsive waveform back-project detection method of lassie (
https://gitext.gfz-potsdam.de/heimann/lassie), and then refines P- and S-wave arrival times 
to produce a robust relocation. For relocation, AutoQuakePycker wraps the NonLinLoc package 
(e.g. Lomax et al., 2009). 

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   modules

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. AutoQuakePycker documentation master file, created by
   sphinx-quickstart on Fri Nov 15 14:49:15 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

AutoQuakePycker
===========================================

AutoQuakePycker is a Python package for the automatic and iterative picking and relocation
of local seismic events. 
AutoQuakePycker takes a first-guess catalogue of earthquake locations (I recommend the very
easy-to-use and comprehnsive waveform back-project detection method of lassie 
(https://gitext.gfz-potsdam.de/heimann/lassie), and then refines P- and S-wave arrival times 
to produce a robust relocation. For relocation, AutoQuakePycker wraps the NonLinLoc package 
(e.g. Lomax et al., 2009). All waveform processing is based on the ObsPy package.

The picking is based on a Kurtosis characteristic function (KCF) which improves pick precision by
computing the KCF over several frequency bandwiths, window sizes and smoothing parameters
(Baillard et al., 2009, BSSA, doi: 10.1785/0120120347).
Bad picks are refined and removed using a clustering procedure (Baillard et al., 2014), as well
as using iterative outlier rejection and distance-dependent residual goals (Sippl et al., 2013, JGR). 
Please see the section below for a full description of the implemented workflow (some aspects
are still work in progress!).

AutoQuakePycker splits the input catalog into chunks to be run as different processes on
multi CPUs for efficient processing of large datasets.

Output hypocentres, arrival times, picks, magnitudes, and polarities
are given in a QUAKEML file per event.

Pre-requisites
--------------
Python functions (all available using Anconda - recommended method):

* NumPy

* SciPy

* matplotlib

* ObsPy (http://www.obspy.org)

Other software:

* NonLinLoc (available to download from http://alomax.free.fr/nlloc/)

In your run directory, you will need to provide a directory called "NLLOC_run" which contains
the (pre-made) travel-time grids and control file comprising the nessecary statements for relocation.
Time2EQ statements are also needed for computing the predicted travel-times for refining the picks.

Detailed description of workflow
----------------

#. Read in preliminary event catalogue (current works for a QUAKEML formatted file).

#. Split catalogue into chunks and split onto multiple cores for efficient processing.


.. toctree::
   :maxdepth: 2
   :caption: API docs:

   modules

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

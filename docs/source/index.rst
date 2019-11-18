.. AutoQuakePycker documentation master file, created by
   sphinx-quickstart on Fri Nov 15 14:49:15 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

AutoQuakePycker
===========================================

AutoQuakePycker is a Python package for the automatic and iterative picking, relocation and
computation of focal mechanism [TODO] for local seismic events. 
AutoQuakePycker takes a first-guess catalogue of earthquake locations (I recommend the very
easy-to-use and comprehnsive waveform back-project detection method of lassie 
(https://gitext.gfz-potsdam.de/heimann/lassie), and then refines P- and S-wave arrival times 
to produce a robust relocation. For relocation, AutoQuakePycker wraps the NonLinLoc package 
(e.g. Lomax et al., 2009). All waveform processing is based on the ObsPy package.

The picking is based on a Kurtosis characteristic function (KCF) which improves pick precision by
computing the KCF over several frequency bandwiths, window sizes and smoothing parameters
(Baillard et al., 2009, BSSA, doi: 10.1785/0120120347). This approach results in a greater number 
of, and higher precision, picked arrival times compared to traditional STA/LTA methods.
Bad picks are refined and removed using a clustering procedure (Baillard et al., 2014), as well
as using iterative outlier rejection and distance-dependent residual goals (Sippl et al., 2013, JGR). 
Please see the section below for a full description of the implemented workflow (some aspects
are still work-in-progress!).

Optionally, S-wave picking is carried out on the horizontal component traces rotated to the 
transverse component, and integrated to displacement.

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

#. Split catalogue into chunks and split onto multiple cores for efficient processing, and process
   each event. Make run directory for each process.

#. Read in available seismic waveform data cut out around initial hypocentre guess. (TODO: Allow the option of 
   reading in continuous data from an archive (e.g. SeisComp Data Structure - SDS).

#. Get theoretical travel times based on initial hypocentre guess and velocity model for stations which have data.

#. For each P and S phase, trim data in given window around theoretical arrival time. Length of window depends 
   on how confident you can be on the initial location.

#. For each filter bandpass and Kurtosis window length, compute the KCF. 

.. toctree::
   :maxdepth: 2
   :caption: API docs:

   modules

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

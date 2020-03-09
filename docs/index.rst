RAPIDpy
=======

RAPIDpy is a python interface for RAPID that assists to prepare inputs,
runs the RAPID program, and provides post-processing utilities. More
information about installation and the input parameters for RAPID can be
found at http://rapid-hub.org. The source code for RAPID is located at
https://github.com/c-h-david/rapid.

|DOI|

|PyPI version|

|Build Status|

|Coverage Status|

|License (3-Clause BSD)|

.. |DOI| image:: https://zenodo.org/badge/19918/erdc/RAPIDpy.svg
   :target: https://zenodo.org/badge/latestdoi/19918/erdc/RAPIDpy
.. |PyPI version| image:: https://badge.fury.io/py/RAPIDpy.svg
   :target: https://badge.fury.io/py/RAPIDpy
.. |Build Status| image:: https://travis-ci.org/erdc/RAPIDpy.svg?branch=master
   :target: https://travis-ci.org/erdc/RAPIDpy
.. |Coverage Status| image:: https://coveralls.io/repos/github/erdc/RAPIDpy/badge.svg?branch=master
   :target: https://coveralls.io/github/erdc/RAPIDpy
.. |License (3-Clause BSD)| image:: https://img.shields.io/badge/license-BSD%203--Clause-yellow.svg
   :target: https://github.com/erdc-cm/RAPIDpy/blob/master/LICENSE

Contents:

.. toctree::
   :maxdepth: 1

   installation
   gis_stream_network    
   gis_tools
   running_rapid
   lsm_inflow
   rapid_dataset
   postprocessing
   rapid_to_gssha


Credit
------

- The development of RAPIDpy was funded by ERDC.
- The pre-processing GIS tools and the inflow tools in RAPIDpy are based on the ESRI RAPID Toolbox (https://github.com/Esri/python-toolbox-for-rapid).


How the Inflow and GIS tools work:
----------------------------------

Snow, Alan D., Scott D. Christensen, Nathan R. Swain, E. James Nelson,
Daniel P. Ames, Norman L. Jones, Deng Ding, Nawajish S. Noman, Cedric H.
David, Florian Pappenberger, and Ervin Zsoter, 2016. A High-Resolution
National-Scale Hydrologic Forecast System from a Global Ensemble Land
Surface Model. *Journal of the American Water Resources Association
(JAWRA)* 1-15, DOI: 10.1111/1752-1688.12434
https://onlinelibrary.wiley.com/doi/full/10.1111/1752-1688.12434

Snow, Alan Dee, "A New Global Forecasting Model to Produce
High-Resolution Stream Forecasts" (2015). All Theses and Dissertations.
Paper 5272.
http://scholarsarchive.byu.edu/etd/5272


Publications using RAPIDpy and RAPID:
-------------------------------------

Tavakoly, A. A., A. D. Snow, C. H. David, M. L. Follum, D. R. Maidment, and Z.-L. Yang, (2016)
"Continental-Scale River Flow Modeling of the Mississippi River Basin Using High-Resolution 
NHDPlus Dataset", Journal of the American Water Resources Association (JAWRA) 1-22. 
DOI: 10.1111/1752-1688.12456


Datasets produced using RAPIDpy and RAPID:
------------------------------------------

Ahmad A Tavakoly. (2017). RAPID input files corresponding to the Mississippi River Basin using 
the NHDPlus v2 Dataset [Data set]. Zenodo. http://doi.org/10.5281/zenodo.322886 |Mississippi Dataset DOI|

.. |Mississippi Dataset DOI| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.322886.svg
   :target: https://doi.org/10.5281/zenodo.322886


Other tools to prepare input for RAPID
---------------------------------------

- For ESRI users: https://github.com/Esri/python-toolbox-for-rapid
- Modified version of the ESRI RAPID Toolbox: https://github.com/erdc-cm/python-toolbox-for-rapid
- For the NHDPlus dataset: https://github.com/c-h-david/RRR



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


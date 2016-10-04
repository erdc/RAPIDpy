Postprocessing
==============

Merge RAPID Output
------------------

.. autoclass:: RAPIDpy.postprocess.ConvertRAPIDOutputToCF

Generate qinit from past qout
-----------------------------

RAPIDpy also creates a qinit file from a RAPID qout file. This example
shows how.

.. automethod:: RAPIDpy.rapid.RAPID.generate_qinit_from_past_qout
    :noindex:                                                

Generate seasonal qinit from past qout
--------------------------------------

.. automethod:: RAPIDpy.rapid.RAPID.generate_seasonal_intitialization
    :noindex:                                                

Goodness of Fit
---------------

To check how well your simulation performed versus observations, these
functions can help you.

.. autofunction:: RAPIDpy.postprocess.find_goodness_of_fit_csv

.. autofunction:: RAPIDpy.postprocess.find_goodness_of_fit


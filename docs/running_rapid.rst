Running RAPID
=============

Tutorial
--------

Step 1: Initialize the RAPID manager class.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
-  First, add the path to you rapid executable location.
-  Next, you need to either tell it how many processors to use using the
   *num\_processors* input variable or to use all available processors
   set *use\_all\_processors* to true.
-  After that, add any other parameters you would like to use that would
   normally be in the rapid namelist file (this is case sensitive).

.. autoclass:: RAPIDpy.rapid.RAPID


Step 2 (optional): Add/update additional namelist parameters later
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automethod:: RAPIDpy.rapid.RAPID.update_parameters

Step 3 (optional): Update reach number data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automethod:: RAPIDpy.rapid.RAPID.update_reach_number_data

Step 4 (optional): Update simulation runtime data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automethod:: RAPIDpy.rapid.RAPID.update_simulation_runtime

Step 5: Run RAPID
~~~~~~~~~~~~~~~~~

.. automethod:: RAPIDpy.rapid.RAPID.run

Step 6 (optional): Convert RAPID output to be CF Compliant
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automethod:: RAPIDpy.rapid.RAPID.make_output_CF_compliant


Full API Description
--------------------

.. autoclass:: RAPIDpy.rapid.RAPID
    :members:
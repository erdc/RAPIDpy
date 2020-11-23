Installation
============

Step 1: Install RAPID
---------------------

Before Installation Steps:
~~~~~~~~~~~~~~~~~~~~~~~~~~

Ubuntu:
^^^^^^^

::

    $ sudo apt-get install gfortran g++

RedHat/CentOS:
^^^^^^^^^^^^^^

::

    $ sudo yum install gcc-c++ gcc-gfortran

Windows with Cygwin:
^^^^^^^^^^^^^^^^^^^^

Downloaded Cygwin (64-bit) (https://www.cygwin.com/) with these
dependencies:

- gcc-core
- gcc-fortran
- gcc-g++
- gdb
- git
- make
- time
- wget

Installation Steps:
~~~~~~~~~~~~~~~~~~~

Manual:
^^^^^^^

-  See: http://rapid-hub.org

Bash:
^^^^^

1. Clone RAPID repository::

    $ git clone https://github.com/c-h-david/rapid.git

2. Install Prereqs::

    $ cd rapid
    $ chmod u+x rapid_install_prereqs.sh
    $ ./rapid_install_prereqs.sh

3. Append *source rapid_specify_varpath.sh* to the ~/.bashrc or ~/.bash_profile::

    source /path/to/cloned/rapid/rapid_specify_varpath.sh

4. Restart Terminal

5. Build RAPID::

    $ cd rapid/src
    $ make rapid

Step 2: Install RAPIDpy
-----------------------

Due to the dependencies required, we recommend using Anaconda or Miniconda.
They can be downloaded from https://www.continuum.io/downloads
or from https://conda.io/miniconda.html.

~~~~~~~~~~~~~~~~~~~~~~

This is how you get the most up-to-date version of the code.

.. note:: If you don't have git, you can download the code from https://github.com/erdc-cm/RAPIDpy

::

    $ git clone https://github.com/erdc-cm/RAPIDpy.git
    $ cd RAPIDpy
    $ conda env create -f rapidpy_env.yml
    $ conda activate rapidpy_env
    $ python setup.py install

To develop on the latest version:

::

    $ git clone https://github.com/erdc-cm/RAPIDpy.git
    $ cd RAPIDpy
    $ conda env create -f rapidpy_env.yml
    $ conda activate rapidpy_env
    $ python setup.py develop

Installation
============

Step 1: Install RAPID
---------------------

Before Installation Steps:
~~~~~~~~~~~~~~~~~~~~~~~~~~

Ubuntu:
^^^^^^^

::

    $ sudo apt-get install gfortran g++ openmpi-bin

RedHat/CentOS:
^^^^^^^^^^^^^^

::

    $ sudo yum install gcc-c++ gcc-gfortran openmpi

Windows with Cygwin:
^^^^^^^^^^^^^^^^^^^^

Downloaded Cygwin (64-bit) (https://www.cygwin.com/) with these
dependencies:

- dos2unix 
- gcc-core 
- gcc-fortran 
- gcc-g++ 
- gdb 
- git
- make 
- netcdf 
- openmpi 
- time 
- wget

Installation Steps:
~~~~~~~~~~~~~~~~~~~

Manual:
^^^^^^^

-  See: http://rapid-hub.org

Bash:
^^^^^

1. Install Prereqs::
    
    $ wget https://raw.githubusercontent.com/snowman2/rapid/master/rapid_install_prereqs.sh
    $ chmod u+x rapid_install_prereqs.sh
    $ ./rapid_install_prereqs.sh

2. Append this to the ~/.bashrc or ~/.bash_profile::

    INSTALLZ_DIR=$HOME/installz
    #-------------------------------------------------------------------------------
    #Exporting environment variables
    #-------------------------------------------------------------------------------
    export TACC_NETCDF_LIB=$INSTALLZ_DIR/netcdf-3.6.3-install/lib
    export TACC_NETCDF_INC=$INSTALLZ_DIR/netcdf-3.6.3-install/include
    export PETSC_DIR=$INSTALLZ_DIR/petsc-3.6.2
    export PETSC_ARCH='linux-gcc-c'
    
    #-------------------------------------------------------------------------------
    #Exporting directories with library-related executables to $PATH
    #-------------------------------------------------------------------------------
    export PATH=$PATH:$PETSC_DIR/$PETSC_ARCH/bin
    export PATH=$PATH:$INSTALLZ_DIR/netcdf-3.6.3-install/bin

3. Restart Terminal

4. Clone RAPID repository and make::

    $ git clone https://github.com/c-h-david/rapid.git
    $ cd rapid/src
    $ make rapid 

Step 2: Install Python Packages
-------------------------------

Method 1: Use Anaconda
~~~~~~~~~~~~~~~~~~~~~~

Linux/Mac
^^^^^^^^^

See: https://github.com/erdc-cm/RAPIDpy/blob/master/.travis.yml

Download Miniconda
''''''''''''''''''

Linux
     

::

    $ wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O miniconda.sh

Mac
   

::

    $ curl -o miniconda.sh https://repo.continuum.io/miniconda/Miniconda2-latest-MacOSX-x86_64.sh

Install Python packages with all dependencies
'''''''''''''''''''''''''''''''''''''''''''''

::

    $ chmod +x miniconda.sh
    $ ./miniconda.sh -b
    $ export PATH=$HOME/miniconda2/bin:$PATH
    $ conda update --yes conda python
    $ conda create --name rapid python=2
    $ source activate rapid
    $ conda install --yes nose numpy scipy netCDF4 gdal shapely pyproj
    $ conda install --yes -c conda-forge rtree
    $ source deactivate rapid
    $ source activate rapid

Windows
^^^^^^^

Download & Install Miniconda
''''''''''''''''''''''''''''

-  Go to: http://conda.pydata.org/miniconda.html
-  Download and run Windows Python 2 version installer
-  Install at
   C:\\Users\\YOUR_USERNAME\\Miniconda2
   or wherever you want
-  Make default python and export to path

Install all dependencies
''''''''''''''''''''''''

Open CMD terminal:

::

    > conda update --yes conda python
    > conda create --name rapid python=2
    > activate rapid
    > conda install --yes nose numpy scipy netCDF4 gdal pyproj pytz python-dateutil
    > conda install --yes -c conda-forge rtree
    > conda install --yes -c scitools shapely
    > deactivate 
    > activate rapid

Method 2: Manual install
~~~~~~~~~~~~~~~~~~~~~~~~

2a: Install netCDF4
^^^^^^^^^^^^^^^^^^^

On Ubuntu:
''''''''''

::

    $ sudo apt-get install python-dev zlib1g-dev libhdf5-serial-dev libnetcdf-dev

On Redhat/CentOS 7:
'''''''''''''''''''

::

    $ sudo yum install netcdf4-python python-devel hdf5-devel netcdf-devel

If you are on RHEL 7 and having troubles, add the epel repo:

::

    $ wget https://dl.fedoraproject.org/pub/epel/epel-release-latest-7.noarch.rpm
    $ sudo rpm -Uvh epel-release-7*.rpm

If you are on CentOS 7 and having troubles, add the epel repo:

::

    $ sudo yum install epel-release

Then install packages listed above.

On OSX:
'''''''

::

    $ brew install homebrew/science/netcdf

2b: (Optional) Install GIS Dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you want to use the GIS preprocessing tools, this section helps you
install the dependencies.

Install GDAL/GEOS/SCIPY:
''''''''''''''''''''''''

Ubuntu:
       

::

    $ sudo apt-get install gdal-bin libproj-dev libgeos-dev python-scipy

RedHat/CentOS:
              

::

    $ sudo yum install gdal proj-devel geos scipy

Install Rtree:
''''''''''''''

See: http://toblerity.org/rtree/install.html

Install spatial python libraries
''''''''''''''''''''''''''''''''

::

    # pip install shapely pyproj gdal rtree

Step 3: Install RAPIDpy
-----------------------

To get the latest stable version:

::

    $ pip install RAPIDpy

To install the latest version:

.. note:: If you don't have git, you can download the code from https://github.com/erdc-cm/RAPIDpy

::

    $ git clone https://github.com/erdc-cm/RAPIDpy.git
    $ cd RAPIDpy
    $ python setup.py install

To develop on the latest version:

::

    $ git clone https://github.com/erdc-cm/RAPIDpy.git
    $ cd RAPIDpy
    $ python setup.py develop

Note: If installing on system, use:

::

    $ sudo su
    # (install command here)
    # exit

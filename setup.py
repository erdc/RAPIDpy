from setuptools import setup, find_packages

setup(
    name='RAPIDpy',
    version='2.2.2',
    description='Python interface for RAPID (rapid-hub.org)',
    long_description='RAPIDpy is a python interface for RAPID that assists to prepare inputs, runs the RAPID program,'
                     ' and provides post-processing utilities (https://github.com/erdc/RAPIDpy). More information '
                     'about installation and the input parameters for RAPID can be found at http://rapid-hub.org. The'
                     ' source code for RAPID is located at https://github.com/c-h-david/rapid/. \n\n'
                     '.. image:: https://zenodo.org/badge/19918/erdc/RAPIDpy.svg \n'
                     '   :target: https://zenodo.org/badge/latestdoi/19918/erdc/RAPIDpy',
    keywords='RAPID',
    author='Alan Dee Snow',
    author_email='alan.d.snow@usace.army.mil',
    url='https://github.com/erdc/RAPIDpy',
    download_url='https://github.com/erdc/RAPIDpy/tarballs/2.2.2',
    license='BSD 3-Clause',
    packages=find_packages(),
    package_data={'': ['gis/lsm_grids/*.nc']},
    install_requires=['numpy', 'netCDF4', 'python-dateutil', 'pytz', 'requests'],
    classifiers=[
                'Intended Audience :: Developers',
                'Intended Audience :: Science/Research',
                'Operating System :: OS Independent',
                'Programming Language :: Python',
                'Programming Language :: Python :: 2',
                'Programming Language :: Python :: 2.7',
                ],
)

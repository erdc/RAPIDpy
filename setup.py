from setuptools import setup, find_packages

setup(
    name='RAPIDpy',
    version='2.4.3',
    description='Python interface for RAPID (rapid-hub.org)',
    long_description='RAPIDpy is a python interface for RAPID that assists '
                     'to prepare inputs, runs the RAPID program, and provides '
                     'post-processing utilities (http://rapidpy.readthedocs.io).'
                     'More information about installation and the input '
                     'parameters for RAPID can be found at http://rapid-hub.org.'
                     ' The source code for RAPID is located at '
                     'https://github.com/c-h-david/rapid/. \n\n'
                     '.. image:: https://zenodo.org/badge/19918/erdc-cm/RAPIDpy.svg \n'
                     '   :target: https://zenodo.org/badge/latestdoi/19918/erdc-cm/RAPIDpy',
    keywords='RAPID',
    author='Alan Dee Snow',
    author_email='alan.d.snow@usace.army.mil',
    url='https://github.com/erdc-cm/RAPIDpy',
    license='BSD 3-Clause',
    packages=find_packages(),
    package_data={'': ['gis/lsm_grids/*.nc']},
    install_requires=['future', 'numpy', 'netcdf4', 'pandas',
                      'python-dateutil', 'pytz', 'requests'],
    classifiers=[
                'Intended Audience :: Developers',
                'Intended Audience :: Science/Research',
                'Operating System :: OS Independent',
                'Programming Language :: Python',
                'Programming Language :: Python :: 2',
                'Programming Language :: Python :: 2.7',
                'Programming Language :: Python :: 3',
                'Programming Language :: Python :: 3.5',
                ],
    extras_require={
        'tests': [
            'pytest',
            'pytest-cov',
        ],
        'docs': [
            'mock',
            'sphinx',
            'sphinx_rtd_theme',
            'sphinxcontrib-napoleon',
        ]
      },
)

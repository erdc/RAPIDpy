from setuptools import setup

setup(
    name='RAPIDpy',
    version='2.1.1',
    description='RAPIDpy is a python interface for RAPID that assists to prepare inputs, runs the RAPID program,'
                ' and provides post-processing utilities (https://github.com/erdc-cm/RAPIDpy). More information '
                'about installation and the input parameters for RAPID can be found at http://rapid-hub.org. The'
                ' source code for RAPID is located at https://github.com/c-h-david/rapid/.',
    keywords='RAPID',
    author='Alan Dee Snow',
    author_email='alan.d.snow@usace.army.mil',
    url='https://github.com/erdc-cm/RAPIDpy',
    download_url='https://github.com/erdc-cm/RAPIDpy/tarballs/2.1.1',
    license='BSD-3 Clause',
    packages=['RAPIDpy'],
    install_requires=['python-dateutil', 'netCDF4', 'numpy', 'pytz', 'requests'],
)

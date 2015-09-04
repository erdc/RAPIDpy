from setuptools import setup

setup(
    name='RAPIDpy',
    version='1.0.3',
    description='Python scripting interface for the RAPID progam. More information about installation and the input parameters can be found at http://rapid-hub.org. The source code for RAPID is located at https://github.com/c-h-david/rapid/.',
    keywords='RAPID',
    author='Alan Dee Snow',
    author_email='alan.d.snow@usace.army.mil',
    url='https://github.com/erdc-cm/RAPIDpy',
    download_url='https://github.com/erdc-cm/RAPIDpy/tarballs/1.0.3',
    license='MIT',
    packages=['RAPIDpy'],
    install_requires=['netCDF4', 'numpy', 'requests'],
)

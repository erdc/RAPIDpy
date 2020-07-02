# -*- coding: utf-8 -*-
#
#  __init__.py
#  pangaea
#
#  Created by Alan D Snow, 2017.
#  BSD 3-Clause
"""pangaea
    Module for reading in land surface model data with xarray.
"""
from .xlsm import LSMGridReader
from .read import open_mfdataset
from .log import log_to_console, log_to_file
from .meta import version

__version__ = version()

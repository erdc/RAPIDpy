#
# Originally from quest python library
# Adapted for pangaea
#
# License: BSD 3-Clause

"""pangaea.log
This module is for logging with pangaea.
Documentation can be found at `_pangaea Documentation HOWTO`_.

.. _pangaea Documentation HOWTO:
   https://github.com/snowman2/pangaea
"""
# default modules
import logging
import os
# external modules
import appdirs
# local modules
from .meta import version

LOGGER = logging.getLogger('pangaea')
LOGGER.addHandler(logging.NullHandler())
LOGGER.propagate = False

DEFAULT_LOG_DIR = appdirs.user_log_dir('pangaea', 'logs')
DEFAULT_LOG_FILE = os.path.join(DEFAULT_LOG_DIR, 'pangaea.log')


def log_to_console(status=True, level=None):
    """Log events to  the console.

    Args:
        status (bool, Optional, Default=True)
            whether logging to console should be turned on(True) or off(False)
        level (string, Optional, Default=None) :
            level of logging; whichever level is chosen all higher levels
            will be logged.
            See: https://docs.python.org/2/library/logging.html#levels
      """

    if status:
        if level is not None:
            LOGGER.setLevel(level)

        console_handler = logging.StreamHandler()
        # create formatter
        formatter = logging.Formatter('%(levelname)s-%(name)s: %(message)s')
        # add formatter to handler
        console_handler.setFormatter(formatter)
        LOGGER.addHandler(console_handler)

        LOGGER.info("pangaea %s", version())

    else:
        for handle in LOGGER.handlers:
            if type(handle).__name__ == 'StreamHandler':
                LOGGER.removeHandler(handle)


def log_to_file(status=True, filename=DEFAULT_LOG_FILE, level=None):
    """Log events to a file.

    Args:
        status (bool, Optional, Default=True)
            whether logging to file should be turned on(True) or off(False)
        filename (string, Optional, Default=None) :
            path of file to log to
        level (string, Optional, Default=None) :
            level of logging; whichever level is chosen all higher levels
            will be logged.
            See: https://docs.python.org/2/library/logging.html#levels
      """

    if status:
        if level is not None:
            LOGGER.setLevel(level)

        try:
            os.mkdir(os.path.dirname(filename))
        except OSError:
            pass

        file_handler = logging.FileHandler(filename)
        # create formatter
        fomat_str = '%(levelname)s-%(name)s: %(message)s'
        formatter = logging.Formatter(fomat_str)
        # add formatter to handler
        file_handler.setFormatter(formatter)
        LOGGER.addHandler(file_handler)

        LOGGER.info("pangaea %s", version())

    else:
        for handle in LOGGER.handlers:
            if type(handle).__name__ == 'FileHandler':
                LOGGER.removeHandler(handle)

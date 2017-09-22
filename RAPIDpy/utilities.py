# -*- coding: utf-8 -*-
"""
   utilities.py
   RAPIDpy

   Created by Alan D. Snow, 2016.
   License BSD-3-Clause
"""
import os
import re

from past.builtins import xrange  # pylint: disable=redefined-builtin


# -----------------------------------------------------------------------------
# HELPER FUNCTIONS
# -----------------------------------------------------------------------------
def case_insensitive_file_search(directory, pattern):
    """
    Looks for file with pattern with case insensitive search
    """
    try:
        return os.path.join(
            directory,
            [filename for filename in os.listdir(directory)
             if re.search(pattern, filename, re.IGNORECASE)][0])
    except IndexError:
        print("{0} not found".format(pattern))
        raise


def partition(lst, n):
    """
        Divide list into n equal parts
    """
    q, r = divmod(len(lst), n)
    indices = [q*i + min(i, r) for i in xrange(n+1)]
    return [lst[indices[i]:indices[i+1]] for i in xrange(n)], \
           [list(xrange(indices[i], indices[i+1])) for i in xrange(n)]


def get_valid_directory_list(input_directory):
    """
    Get a list of folders
    """
    valid_input_directories = []
    for directory in os.listdir(input_directory):
        if os.path.isdir(os.path.join(input_directory, directory)):
            valid_input_directories.append(directory)
        else:
            print("{0} not a directory. Skipping ...".format(directory))
    return valid_input_directories

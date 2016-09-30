# -*- coding: utf-8 -*-
##
##  utilities.py
##  RAPIDpy
##
##  Created by Alan D. Snow.
##  Copyright © 2016 Alan D Snow. All rights reserved.
##

import os
from past.builtins import xrange
import re

#----------------------------------------------------------------------------------------
# HELPER FUNCTIONS
#----------------------------------------------------------------------------------------
def case_insensitive_file_search(directory, pattern):
    """
    Looks for file with pattern with case insensitive search
    """
    try:
        return os.path.join(directory,
                            [filename for filename in os.listdir(directory) \
                             if re.search(pattern, filename, re.IGNORECASE)][0])
    except IndexError:
        print("{0} not found".format(pattern))
        raise

def partition(lst, n):
    """
        Divide list into n equal parts
    """
    q, r = divmod(len(lst), n)
    indices = [q*i + min(i,r) for i in xrange(n+1)]
    return [lst[indices[i]:indices[i+1]] for i in xrange(n)], \
           [list(xrange(indices[i],indices[i+1])) for i in xrange(n)]

def get_valid_watershed_list(input_directory):
    """
    Get a list of folders formatted correctly for watershed-subbasin
    """
    valid_input_directories = []
    for directory in os.listdir(input_directory):
        if os.path.isdir(os.path.join(input_directory, directory)) \
            and len(directory.split("-")) == 2:
            valid_input_directories.append(directory)
        else:
            print("{0} incorrectly formatted. Skipping ...".format(directory))
    return valid_input_directories

def get_watershed_subbasin_from_folder(folder_name):
    """
    Get's the watershed & subbasin name from folder
    """
    input_folder_split = folder_name.split("-")
    watershed = input_folder_split[0].lower()
    subbasin = input_folder_split[1].lower()
    return watershed, subbasin

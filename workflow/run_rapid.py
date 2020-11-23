#! /usr/bin/env python

from RAPIDpy import RAPID
import sys

if __name__=='__main__':
    rapid_namelist_file = sys.argv[1]
    try:
	qinit_file = sys.argv[2]
	qinit = True
    except:
	qinit_file = ''
	qinit = False
 
    rapid_manager = RAPID(
        rapid_executable_location='/home/mgeheran/rapid/rapid/run/rapid',
        use_all_processors=True, Qinit_file=qinit_file, BS_opt_Qinit=qinit)

    rapid_manager.run(rapid_namelist_file=rapid_namelist_file)
        
    print("completed.")



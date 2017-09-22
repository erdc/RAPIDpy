# -*- coding: utf-8 -*-
"""
    rapid.py
    RAPIDpy

    Created by Alan D Snow, 2015.
    License: BSD-3-Clause
"""
from calendar import isleap
from csv import writer as csvwriter
import datetime
from multiprocessing import cpu_count
import os
from subprocess import Popen, PIPE
from time import gmtime

from dateutil.parser import parse
import numpy as np
from requests import get
import xarray

from .dataset import RAPIDDataset
from .helper_functions import csv_to_list, log, open_csv
from .postprocess import ConvertRAPIDOutputToCF


# -----------------------------------------------------------------------------
# Main RAPID Manager Class
# -----------------------------------------------------------------------------
class RAPID(object):
    """
    This class is designed to prepare the rapid_namelist file and run
    the RAPID program. There are also other utilities added.

    Attributes
    ----------
    rapid_executable_location: str, optional
        Path to the RAPID executable location.
    num_processors: int, optional
        Number of procesors to use. Default is 1.
        Overridden if *use_all_processors* is True.
    use_all_processors: bool, optional
        If set to True, the RAPID program will use all available processors.
        Default is False.
    cygwin_bin_location: str, optional
        If using Windows, this is the path to the Cygwin 'bin' directory.
    mpiexec_command: str, optional
        This is the mpi execute commmand. Default is "mpiexec".
    ksp_type: str, optional
        This is the solver type. Default is "richardson".
    **kwargs: str, optional
        Keyword arguments matching the input parameters in the RAPID namelist.


    Linux Example:

    .. code:: python

        from RAPIDpy import RAPID

        rapid_manager = RAPID(
            rapid_executable_location='~/work/rapid/run/rapid'
            use_all_processors=True,
            ZS_TauR=24 * 3600,
            ZS_dtR=15 * 60,
            ZS_TauM=365 * 24 * 3600,
            ZS_dtM=24 * 3600
        )


    Windows with Cygwin Example:

    .. code:: python

        from RAPIDpy import RAPID

        cygwin_exe = 'C:/cygwin64/home/username/work/rapid/run/rapid'
        rapid_manager = RAPID(
            rapid_executable_location=cygwin_exe,
            cygwin_bin_location='C:/cygwin64/bin',
            use_all_processors=True,
            ZS_TauR=24 * 3600,
            ZS_dtR=15 * 60,
            ZS_TauM=365 * 24 * 3600,
            ZS_dtM=24 * 3600
        )

    """
    # pylint: disable=too-many-instance-attributes
    def __init__(self,
                 rapid_executable_location="",
                 num_processors=1,
                 use_all_processors=False,
                 cygwin_bin_location="",
                 mpiexec_command="mpiexec",
                 ksp_type="richardson",
                 **kwargs):
        """
        Initialize the class with variables given by the user
        """
        if os.name == "nt" and not \
                (cygwin_bin_location or os.path.exists(cygwin_bin_location))\
                and rapid_executable_location:
            raise Exception("Required to have cygwin_bin_location set "
                            "if using windows!")

        self._rapid_executable_location = rapid_executable_location
        self._cygwin_bin_location = cygwin_bin_location
        self._cygwin_bash_exe_location = \
            os.path.join(cygwin_bin_location, "bash.exe")
        self._mpiexec_command = mpiexec_command
        self._ksp_type = ksp_type

        # use all processors makes precedent over num_processors arg
        if use_all_processors is True:
            self._num_processors = cpu_count()
        elif num_processors > cpu_count():
            log("Num processors requested exceeded max. Set to max ...",
                "WARNING")
            self._num_processors = cpu_count()
        else:
            self._num_processors = num_processors

        # ---------------------------------------------------------------------
        # Runtime options
        # ---------------------------------------------------------------------
        self.BS_opt_Qinit = False
        # !.false. --> no read initial flow    .true. --> read initial flow
        self.BS_opt_Qfinal = False
        # !.false. --> no write final flow     .true. --> write final flow
        self.BS_opt_dam = False
        # !.false. --> no dam model used       .true. --> dam model used
        self.BS_opt_for = False
        # !.false. --> no forcing              .true. --> forcing
        self.BS_opt_influence = False
        # !.false. --> no output influence     .true. --> output influence
        self.IS_opt_routing = 1
        # !1      --> matrix-based Muskingum  2      --> traditional Muskingum
        # !3      --> Transbnd. matrix-based
        self.IS_opt_run = 1
        # !1      --> regular run             2      --> parameter optimization
        self.IS_opt_phi = 1
        # !1      --> phi1                    2      --> phi2
        # ---------------------------------------------------------------------
        # Temporal information
        # ---------------------------------------------------------------------
        # NOTE: ALL TIME IN SECONDS!
        # ALWAYS USED
        self.ZS_TauR = 0
        # duration of routing procedure (time step of runoff data)
        self.ZS_dtR = 0
        # internal routing time step
        # ONLY FOR REGULAR RUN
        self.ZS_TauM = 0
        # total simulation time
        self.ZS_dtM = 0
        # input time step
        # ONLY FOR OPTIMIZATION RUN
        self.ZS_TauO = 0
        # total optimization time
        self.ZS_dtO = 0
        # observation time step
        # FORCING MODE (replace some values with observations)
        self.ZS_dtF = 0
        # time step of forcing data
        # ---------------------------------------------------------------------
        # Domain in which input data is available
        # ---------------------------------------------------------------------
        self.IS_riv_tot = 0
        # number of river reaches in rapid connect file
        self.rapid_connect_file = ''
        # path to rapid_connect file
        self.IS_max_up = 0
        # maximum number of ustream segments
        self.Vlat_file = ''
        # path to runoff file
        # ---------------------------------------------------------------------
        # Domain in which model runs
        # ---------------------------------------------------------------------
        self.IS_riv_bas = 0
        # number of river reaches in subbasin
        self.riv_bas_id_file = ''
        # subbasin reach id file
        # ---------------------------------------------------------------------
        # Initial instantaneous flow file
        # ---------------------------------------------------------------------
        self.Qinit_file = ''
        # initial flow file (same order as rapid_connect)
        # ---------------------------------------------------------------------
        # Final instantaneous flow file
        # ---------------------------------------------------------------------
        self.Qfinal_file = ''
        # path to output final flow file
        # ---------------------------------------------------------------------
        # Available dam data
        # ---------------------------------------------------------------------
        self.IS_dam_tot = 0
        # number of dams
        self.dam_tot_id_file = ''
        # ids of dam location
        # ---------------------------------------------------------------------
        # Dam data used
        # ---------------------------------------------------------------------
        self.IS_dam_use = 0
        # number in subset of dam data to use
        self.dam_use_id_file = ''
        # ids of subset of dams
        # ---------------------------------------------------------------------
        # Available forcing data
        # ---------------------------------------------------------------------
        self.IS_for_tot = 0
        self.for_tot_id_file = ''
        self.Qfor_file = ''
        # ---------------------------------------------------------------------
        # Forcing data used as model runs
        # ---------------------------------------------------------------------
        self.IS_for_use = 0
        self.for_use_id_file = ''
        # ---------------------------------------------------------------------
        # File where max (min) of absolute values of b (QoutR) are stored
        # ---------------------------------------------------------------------
        self.babsmax_file = ''
        self.QoutRabsmin_file = ''
        self.QoutRabsmax_file = ''
        # ---------------------------------------------------------------------
        # Regular model run
        # ---------------------------------------------------------------------
        self.k_file = ''
        self.x_file = ''
        self.Qout_file = ''
        # ---------------------------------------------------------------------
        # Optimization
        # ---------------------------------------------------------------------
        self.ZS_phifac = 0
        # ---------------------------------------------------------------------
        # Routing parameters
        # ---------------------------------------------------------------------
        self.kfac_file = ''
        self.xfac_file = ''
        self.ZS_knorm_init = 0
        self.ZS_xnorm_init = 0
        # ---------------------------------------------------------------------
        # Gage observations
        # ---------------------------------------------------------------------
        self.IS_obs_tot = 0
        self.obs_tot_id_file = ''
        self.Qobs_file = ''
        self.Qobsbarrec_file = ''
        self.IS_obs_use = 0
        self.obs_use_id_file = ''
        self.IS_strt_opt = 0

        self.update_parameters(**kwargs)

    def _get_cygwin_path(self, windows_path):
        """
        Convert windows path to cygpath
        """
        conv_cmd = [os.path.join(self._cygwin_bin_location, "cygpath.exe"),
                    "-u", windows_path]
        process = Popen(conv_cmd,
                        stdout=PIPE, stderr=PIPE, shell=False)
        out, err = process.communicate()
        if err:
            print(err)
            raise Exception(err)

        return out.strip()

    def _create_symlink_cygwin(self, initial_path, final_path):
        """
        Use cygqin to generate symbolic link
        """
        symlink_cmd = [os.path.join(self._cygwin_bin_location, "ln.exe"),
                       "-s", self._get_cygwin_path(initial_path),
                       self._get_cygwin_path(final_path)]
        process = Popen(symlink_cmd,
                        stdout=PIPE, stderr=PIPE, shell=False)
        out, err = process.communicate()
        if err:
            print(err)
            raise Exception(err)

        return out.strip()

    def _dos2unix_cygwin(self, file_path):
        """
        Use cygwin to convert file to unix format
        """
        dos2unix_cmd = \
            [os.path.join(self._cygwin_bin_location, "dos2unix.exe"),
             self._get_cygwin_path(file_path)]
        process = Popen(dos2unix_cmd,
                        stdout=PIPE, stderr=PIPE, shell=False)
        process.communicate()

    def update_parameters(self, **kwargs):
        """
        You can add or update rapid namelist parameters by using the name of
        the variable in the rapid namelist file (this is case sensitive).

        Parameters
        ----------
        **kwargs: str, optional
            Keyword arguments matching the input parameters
            in the RAPID namelist.


        Example:

        .. code:: python

            from RAPIDpy import RAPID

            rapid_manager = RAPID(
                rapid_executable_location='~/work/rapid/run/rapid'
                use_all_processors=True,
                ZS_TauR=24 * 3600,
                ZS_dtR=15 * 60,
                ZS_TauM=365 * 24 * 3600,
                ZS_dtM=24 * 3600
            )

            rapid_manager.update_parameters(
                rapid_connect_file='../rapid-io/input/rapid_connect.csv',
                Vlat_file='../rapid-io/input/m3_riv.nc',
                riv_bas_id_file='../rapid-io/input/riv_bas_id.csv',
                k_file='../rapid-io/input/k.csv',
                x_file='../rapid-io/input/x.csv',
                Qout_file='../rapid-io/output/Qout.nc',
            )

        """
        # set arguments based off of user input
        for key, value in list(kwargs.items()):
            if key in dir(self) and not key.startswith('_'):
                setattr(self, key, value)
            else:
                log("Invalid RAPID parameter %s." % key,
                    "ERROR")

    def update_reach_number_data(self):
        """
        Update the reach number data for the namelist based on input files.

        .. warning:: You need to make sure you set *rapid_connect_file*
                     and *riv_bas_id_file* before running this function.


        Example:

        .. code:: python

            from RAPIDpy import RAPID

            rapid_manager = RAPID(
                rapid_connect_file='../rapid-io/input/rapid_connect.csv',
                riv_bas_id_file='../rapid-io/input/riv_bas_id.csv',
            )

            rapid_manager.update_reach_number_data()


        Example with forcing data:

        .. code:: python

            from RAPIDpy import RAPID

            rapid_manager = RAPID(
                rapid_connect_file='../rapid-io/input/rapid_connect.csv',
                riv_bas_id_file='../rapid-io/input/riv_bas_id.csv',
                Qfor_file='../rapid-io/input/qfor_file.csv',
                for_tot_id_file='../rapid-io/input/for_tot_id_file.csv',
                for_use_id_file='../rapid-io/input/for_use_id_file.csv',
                ZS_dtF=3*60*60,
                BS_opt_for=True
            )

            rapid_manager.update_reach_number_data()

        """
        if not self.rapid_connect_file:
            log("Missing rapid_connect_file. "
                "Please set before running this function ...",
                "ERROR")

        if not self.riv_bas_id_file:
            log("Missing riv_bas_id_file. "
                "Please set before running this function ...",
                "ERROR")

        # get rapid connect info
        rapid_connect_table = np.loadtxt(self.rapid_connect_file,
                                         ndmin=2, delimiter=",", dtype=int)

        self.IS_riv_tot = int(rapid_connect_table.shape[0])
        self.IS_max_up = int(rapid_connect_table[:, 2].max())

        # get riv_bas_id info
        riv_bas_id_table = np.loadtxt(self.riv_bas_id_file,
                                      ndmin=1, delimiter=",",
                                      usecols=(0,), dtype=int)
        self.IS_riv_bas = int(riv_bas_id_table.size)

        # add the forcing files
        if not self.for_tot_id_file:
            self.IS_for_tot = 0
            log("Missing for_tot_id_file. Skipping ...",
                "WARNING")
        else:
            # get riv_bas_id info
            for_tot_id_table = np.loadtxt(self.for_tot_id_file,
                                          ndmin=1, delimiter=",",
                                          usecols=(0,), dtype=int)
            self.IS_for_tot = int(for_tot_id_table.size)

        if not self.for_use_id_file:
            self.IS_for_use = 0
            log("Missing for_use_id_file. Skipping ...",
                "WARNING")
        else:
            # get riv_bas_id info
            for_use_id_table = np.loadtxt(self.for_use_id_file,
                                          ndmin=1, delimiter=",",
                                          usecols=(0,), dtype=int)
            self.IS_for_use = int(for_use_id_table.size)

    def update_simulation_runtime(self):
        """
        Updates the total simulation duration from
        the m3 file (Vlat_file) and the time step (ZS_TauR).

        .. warning:: You need to set the m3 file (Vlat_file) and the
                     time step (ZS_TauR) before runnning this function.


        Example:

        .. code:: python

            from RAPIDpy import RAPID

            rapid_manager = RAPID(
                Vlat_file='../rapid-io/input/m3_riv.csv',
                ZS_TauR=3*3600,
            )

            rapid_manager.update_simulation_runtime()
        """
        if not self.Vlat_file or not os.path.exists(self.Vlat_file):
            log("Need Vlat_file to proceed ...",
                "ERROR")

        if self.ZS_TauR <= 0:
            log("Missing routing time step ...",
                "ERROR")

        try:
            self.ZS_TauR = int(self.ZS_TauR)
        except ValueError:
            log("Invalid routing time step: {0} ...".format(self.ZS_TauR),
                "ERROR")

        with RAPIDDataset(self.Vlat_file) as m3_nc:
            self.ZS_TauM = m3_nc.size_time*self.ZS_TauR
            self.ZS_TauO = m3_nc.size_time*self.ZS_TauR

    def generate_namelist_file(self, rapid_namelist_file):
        """
        Generate rapid_namelist file.

        Parameters
        ----------
        rapid_namelist_file: str
            Path of namelist file to generate from
            parameters added to the RAPID manager.
        """
        log("Generating RAPID namelist file ...",
            "INFO")
        try:
            os.remove(rapid_namelist_file)
        except OSError:
            pass

        with open(rapid_namelist_file, 'w') as new_file:
            new_file.write('&NL_namelist\n')
            for attr, value in sorted(list(self.__dict__.items())):
                if not attr.startswith('_'):
                    if attr.startswith('BS'):
                        new_file.write("{0} = .{1}.\n"
                                       .format(attr, str(value).lower()))
                    elif isinstance(value, int):
                        new_file.write("%s = %s\n" % (attr, value))
                    else:
                        if value:
                            if os.name == "nt":
                                # if windows generate file with cygpath
                                value = self._get_cygwin_path(value)
                            new_file.write("%s = \'%s\'\n" % (attr, value))
            new_file.write("/\n")

    def update_namelist_file(self, rapid_namelist_file,
                             new_namelist_file=None):
        """
        Update existing namelist file with new parameters

        Parameters
        ----------
        rapid_namelist_file: str
            Path of namelist file to use in the simulation. It will be
            updated with any parameters added to the RAPID manager.
        new_namelist_file: str, optional
            Path to output the updated namelist file.
        """
        if os.path.exists(rapid_namelist_file) and rapid_namelist_file:
            log("Adding missing inputs from RAPID input file ...",
                "INFO")
            with open(rapid_namelist_file, 'r') as old_file:
                for line in old_file:
                    line = line.strip()
                    if not line[:1].isalpha() or not line:
                        continue
                    line_split = line.split("=")
                    attr = line_split[0].strip()
                    value = None
                    if len(line_split) > 1:
                        value = line_split[1].strip()\
                            .replace("'", "").replace('"', "")
                        # convert integers to integers
                        try:
                            value = int(value)
                        except ValueError:
                            pass
                        # remove dots from beginning & end of value
                        if attr.startswith('BS'):
                            value = value.replace(".", "")
                    # add attribute if exists
                    if attr in dir(self) and not attr.startswith('_'):
                        # set attribute if not set already
                        if not getattr(self, attr):
                            setattr(self, attr, value)
                    else:
                        log("Invalid argument {0}. Skipping ...".format(attr),
                            "INFO")

            if new_namelist_file is None:
                new_namelist_file = rapid_namelist_file

            self.generate_namelist_file(new_namelist_file)
        else:
            log("RAPID namelist file to update not found.",
                "ERROR")

    def make_output_cf_compliant(self,
                                 simulation_start_datetime,
                                 comid_lat_lon_z_file="",
                                 project_name="Normal RAPID project"):
        """
        This function converts the RAPID output to be CF compliant.
        This will require a *comid_lat_lon_z.csv* file
        (See: :func:`~RAPIDpy.gis.centroid.FlowlineToPoint` to
        generate the file).

        .. note:: It prepends time an initial flow to your simulation from the
                  *qinit_file*. If no qinit file is given, an initial value
                  of zero is added.

        .. warning:: This will delete your original Qout file.

        Parameters
        ----------
        simulation_start_datetime: datetime
            Datetime object with the start date of the simulation.
        comid_lat_lon_z_file: str, optional
            Path to the *comid_lat_lon_z.csv* file. If none given,
            spatial information will be skipped.
        project_name: str, optional
            Name of project to add to the RAPID output file.


        Example:

        .. code:: python

            from RAPIDpy import RAPID

            rapid_manager = RAPID(
                rapid_executable_location='~/work/rapid/run/rapid'
                use_all_processors=True,
                ZS_TauR=24*3600,
                ZS_dtR=15*60,
                ZS_TauM=365*24*3600,
                ZS_dtM=24*3600
                rapid_connect_file='../rapid-io/input/rapid_connect.csv',
                Vlat_file='../rapid-io/input/m3_riv.nc',
                riv_bas_id_file='../rapid-io/input/riv_bas_id.csv',
                k_file='../rapid-io/input/k.csv',
                x_file='../rapid-io/input/x.csv',
                Qout_file='../rapid-io/output/Qout.nc'
            )

            rapid_manager.run()

            rapid_manager.make_output_cf_compliant(
                simulation_start_datetime=datetime.datetime(1980, 1, 1),
                comid_lat_lon_z_file='../rapid-io/input/comid_lat_lon_z.csv',
                project_name="ERA Interim Historical flows by US Army ERDC"
            )

        """
        with RAPIDDataset(self.Qout_file) as qout_nc:
            if qout_nc.is_time_variable_valid():
                log("RAPID Qout file already CF compliant ...",
                    "INFO")
                return

        crv = ConvertRAPIDOutputToCF(
            rapid_output_file=self.Qout_file,
            start_datetime=simulation_start_datetime,
            time_step=self.ZS_TauR,
            qinit_file=self.Qinit_file,
            comid_lat_lon_z_file=comid_lat_lon_z_file,
            rapid_connect_file=self.rapid_connect_file,
            project_name=project_name,
            output_id_dim_name='rivid',
            output_flow_var_name='Qout',
            print_debug=False
        )
        crv.convert()

    def run(self, rapid_namelist_file=""):
        """
        Run RAPID program and generate file based on inputs
        This will generate your rapid_namelist file and run RAPID from wherever
        you call this script (your working directory).

        Parameters
        ----------
        rapid_namelist_file: str, optional
            Path of namelist file to use in the simulation.
            It will be updated with any parameters added to the RAPID manager.


        Linux Example:

        .. code:: python

            from RAPIDpy import RAPID

            rapid_manager = RAPID(
                rapid_executable_location='~/work/rapid/src/rapid'
                use_all_processors=True,
            )

            rapid_manager.update_parameters(
                rapid_connect_file='../rapid-io/input/rapid_connect.csv',
                Vlat_file='../rapid-io/input/m3_riv.nc',
                riv_bas_id_file='../rapid-io/input/riv_bas_id.csv',
                k_file='../rapid-io/input/k.csv',
                x_file='../rapid-io/input/x.csv',
                Qout_file='../rapid-io/output/Qout.nc',
            )

            rapid_manager.update_reach_number_data()
            rapid_manager.update_simulation_runtime()
            rapid_manager.run(
                rapid_namelist_file='../rapid-io/input/rapid_namelist')


        Linux Reservoir Forcing Flows Example:

        .. code:: python

            from RAPIDpy import RAPID

            rapid_manager = RAPID(
                rapid_executable_location='~/work/rapid/src/rapid',
                num_processors=4,
                IS_for_tot=4,
                IS_for_use=4,
                for_tot_id_file='../rapid-io/input/dam_id.csv',
                for_use_id_file='../rapid-io/input/dam_id.csv',
                Qfor_file='../rapid-io/input/qout_dams.csv',
                ZS_dtF=86400,
                BS_opt_for=True,
            )

            rapid_manager.run(
                rapid_namelist_file='../rapid-io/input/rapid_namelist_regular')

        Windows with Cygwin Example:

        .. code:: python

            from RAPIDpy import RAPID
            from os import path

            rapid_exe_path = 'C:/cygwin64/home/username/rapid/run/rapid',
            rapid_manager = RAPID(
                rapid_executable_location=rapid_exe_path,
                cygwin_bin_location='C:/cygwin64/bin',
                use_all_processors=True,
                ZS_TauR=24*3600,
                ZS_dtR=15*60,
                ZS_TauM=365*24*3600,
                ZS_dtM=24*3600
            )

            rapid_input = 'C:/cygwin64/home/username/rapid-io/input'
            rapid_output = 'C:/cygwin64/home/username/rapid-io/output'
            rapid_manager.update_parameters(
                rapid_connect_file=path.join(rapid_input, 'rapid_connect.csv'),
                Vlat_file=path.join(rapid_input, 'm3_riv.nc'),
                riv_bas_id_file=path.join(rapid_input, 'riv_bas_id.csv'),
                k_file=path.join(rapid_input, 'k.csv'),
                x_file=path.join(rapid_input, 'x.csv'),
                Qout_file=path.join(rapid_output, 'Qout.nc'),
            )

            rapid_manager.update_reach_number_data()
            rapid_manager.update_simulation_runtime()
            rapid_manager.run()
        """
        if not self._rapid_executable_location:
            log("Missing rapid_executable_location. "
                "Please set before running this function ...",
                "ERROR")

        time_start = datetime.datetime.utcnow()
        temp_rapid_namelist_file = os.path.join(os.getcwd(), "rapid_namelist")

        if not rapid_namelist_file or not os.path.exists(rapid_namelist_file):
            # generate input file if it does not exist
            self.generate_namelist_file(temp_rapid_namelist_file)
        else:
            # update existing file
            self.update_namelist_file(rapid_namelist_file,
                                      temp_rapid_namelist_file)

        local_rapid_executable_location = \
            os.path.join(os.path.dirname(temp_rapid_namelist_file),
                         "rapid_exe_symlink")

        def rapid_cleanup(*args):
            """
            Cleans up the rapid files generated by the process
            """
            for arg in args:
                # remove files
                try:
                    os.remove(arg)
                except OSError:
                    pass

        # create link to RAPID if needed
        temp_link_to_rapid = ""
        # pylint: disable=no-member
        if self._rapid_executable_location != \
                local_rapid_executable_location:
            rapid_cleanup(local_rapid_executable_location)
            if os.name == "nt":
                self._create_symlink_cygwin(self._rapid_executable_location,
                                            local_rapid_executable_location)
            else:
                os.symlink(self._rapid_executable_location,
                           local_rapid_executable_location)
            temp_link_to_rapid = local_rapid_executable_location

        # run RAPID
        log("Running RAPID ...",
            "INFO")
        if os.name == "nt":
            local_rapid_executable_location = \
                self._get_cygwin_path(local_rapid_executable_location)

        # htcondor will not allow mpiexec for single processor jobs
        # this was added for that purpose
        run_rapid_command = [local_rapid_executable_location,
                             "-ksp_type", self._ksp_type]

        if self._num_processors > 1:
            run_rapid_command = [self._mpiexec_command,
                                 "-n", str(self._num_processors)] \
                                + run_rapid_command

        process = Popen(run_rapid_command,
                        stdout=PIPE, stderr=PIPE, shell=False)
        out, err = process.communicate()
        if err:
            rapid_cleanup(temp_link_to_rapid, temp_rapid_namelist_file)
            raise Exception(err)
        else:
            log('RAPID output:',
                "INFO")
            for line in out.split(b'\n'):
                print(line)
        rapid_cleanup(temp_link_to_rapid, temp_rapid_namelist_file)
        log("Time to run RAPID: %s" % (datetime.datetime.utcnow()-time_start),
            "INFO")

    def generate_qinit_from_past_qout(self, qinit_file, time_index=-1,
                                      out_datetime=None):
        """
        Generate qinit from a RAPID qout file

        Parameters
        ----------
        qinit_file: str
            Path to output qinit_file.
        time_index: int, optional
            Index of simulation to generate initial flow file.
            Default is the last index.
        out_datetime: :obj:`datetime.datetime`, optional
            Datetime object containing time of initialization.


        Example:

        .. code:: python

            from RAPIDpy import RAPID

            rapid_manager = RAPID(
                Qout_file='/output_mississippi-nfie/Qout_k2v1_2005to2009.nc',
                rapid_connect_file='/input_mississippi_nfie/rapid_connect.csv'
            )

            rapid_manager.generate_qinit_from_past_qout(
                qinit_file='/input_mississippi_nfie/Qinit_2008_flood.csv',
                time_index=10162
            )

        """
        if not self.Qout_file or not os.path.exists(self.Qout_file):
            log('Missing Qout_file. '
                'Please set before running this function ...',
                "ERROR")

        if not self.rapid_connect_file or not self.rapid_connect_file:
            log('Missing rapid_connect file. '
                'Please set before running this function ...',
                "ERROR")

        log("Generating qinit file from qout file ...",
            "INFO")
        # get information from dataset
        with xarray.open_dataset(self.Qout_file) as qds:
            rivid_array = qds.rivid.values
            if out_datetime is None:
                streamflow_values = qds.isel(time=time_index).Qout.values
            else:
                streamflow_values = qds.sel(time=str(out_datetime)).Qout.values

        log("Reordering data ...",
            "INFO")

        stream_id_array = np.loadtxt(self.rapid_connect_file,
                                     ndmin=1, delimiter=",",
                                     usecols=(0,), dtype=int)
        init_flows_array = np.zeros(stream_id_array.size)
        for riv_bas_index, riv_bas_id in enumerate(rivid_array):
            try:
                data_index = np.where(stream_id_array == riv_bas_id)[0][0]
                init_flows_array[data_index] = streamflow_values[riv_bas_index]
            except IndexError:
                log('riv bas id {0} not found in connectivity list.'
                    .format(riv_bas_id),
                    "WARNING")

        log("Writing to file ...",
            "INFO")
        with open_csv(qinit_file, 'w') as qinit_out:
            for init_flow in init_flows_array:
                qinit_out.write('{0}\n'.format(init_flow))

        self.Qinit_file = qinit_file
        self.BS_opt_Qinit = True
        log("Initialization Complete!",
            "INFO")

    def generate_seasonal_intitialization(
        self,
        qinit_file,
        datetime_start_initialization=datetime.datetime.utcnow()
    ):
        """This creates a seasonal qinit file from a RAPID qout file. This
        requires a simulation Qout file with a longer time period of record and
        to be CF compliant. It takes the average of the current date +- 3 days
        and goes back as far as possible.

        Parameters
        ----------
        qinit_file: str
            Path to output qinit_file.
        datetime_start_initialization: :obj:`datetime.datetime`, optional
            Datetime object with date of simulation to go back through the
            years and get a running average to generate streamflow
            initialization. Default is utcnow.


        Example:

        .. code:: python

            from RAPIDpy.rapid import RAPID

            rapid_manager = RAPID(
                Qout_file='/output_mississippi-nfie/Qout_2000to2015.nc',
                rapid_connect_file='/input_mississippi_nfie/rapid_connect.csv'
            )

            rapid_manager.generate_seasonal_intitialization(
                qinit_file='/input_mississippi_nfie/Qinit_seasonal_avg.csv'
            )
        """
        if not self.Qout_file or not os.path.exists(self.Qout_file):
            log("Missing Qout_file. "
                "Please set before running this function ...",
                "ERROR")

        if not self.rapid_connect_file or not self.rapid_connect_file:
            log("Missing rapid_connect file. "
                "Please set before running this function ...",
                "ERROR")

        day_of_year = datetime_start_initialization.timetuple().tm_yday
        min_day = day_of_year - 3
        max_day = day_of_year + 3

        with RAPIDDataset(self.Qout_file) as qout_hist_nc:
            if not qout_hist_nc.is_time_variable_valid():
                log("File must be CF 1.6 compliant "
                    "with valid time variable ...",
                    "ERROR")

            log("Generating seasonal average qinit file from qout file ...",
                "INFO")

            log("Determining dates with streamflows of interest ...",
                "INFO")

            time_indices = []
            for idx, ttt in enumerate(qout_hist_nc.get_time_array()):
                var_time = gmtime(ttt)
                compare_yday = var_time.tm_yday
                # move day back one past because of leap year adds
                # a day after feb 29 (day 60)
                if isleap(var_time.tm_year) and compare_yday > 60:
                    compare_yday -= 1
                # check if date within range of season
                if min_day <= compare_yday < max_day:
                    time_indices.append(idx)

            if not time_indices:
                log("No time steps found within range ...",
                    "ERROR")

            log("Extracting data ...",
                "INFO")

            streamflow_array = \
                qout_hist_nc.get_qout(time_index_array=time_indices)

            log("Reordering data...",
                "INFO")
            stream_id_array = np.loadtxt(self.rapid_connect_file,
                                         ndmin=1, delimiter=",",
                                         usecols=(0,), dtype=int)
            init_flows_array = np.zeros(stream_id_array.size)
            for riv_bas_index, riv_bas_id in enumerate(
                    qout_hist_nc.get_river_id_array()):
                try:
                    data_index = np.where(stream_id_array == riv_bas_id)[0][0]
                    init_flows_array[data_index] = \
                        np.mean(streamflow_array[riv_bas_index])
                except IndexError:
                    log('riv_bas_id {0} not found in connectivity list.'
                        .format(riv_bas_id),
                        "WARNING")

            log("Writing to file ...",
                "INFO")
            with open_csv(qinit_file, 'w') as qinit_out:
                for init_flow in init_flows_array:
                    qinit_out.write('{}\n'.format(init_flow))

            log("Initialization Complete!",
                "INFO")

    def generate_usgs_avg_daily_flows_opt(self,
                                          reach_id_gage_id_file,
                                          start_datetime,
                                          end_datetime,
                                          out_streamflow_file,
                                          out_stream_id_file):
        """
        Generate daily streamflow file and stream id file required for
        calibration or for substituting flows based on USGS gage ids
        associated with stream ids.

        Parameters
        ----------
        reach_id_gage_id_file: str
            Path to reach_id_gage_id file.
        start_datetime: datetime
            A datetime object with the start date to download data.
        end_datetime: datetime
            A datetime object with the end date to download data.
        out_streamflow_file: str
            The path to output the streamflow file for RAPID.
        out_stream_id_file: str
            The path to output the stream ID file associated with the
            streamflow file for RAPID.


        Example *reach_id_gage_id_file*::

            COMID, USGS_GAGE_ID
            2000, 503944
            ...

        .. warning:: Overuse will get you blocked from downloading data from
                     USGS.

        .. warning:: This code does not clean the data in any way. Thus, you
                     are likely to run into issues if you simply use the raw
                     data.

        .. warning:: The code skips gages that do not have data
                     for the entire time period.


        Simple Example:

        .. code:: python

            import datetime
            from os.path import join
            from RAPIDpy import RAPID

            main_path = "/home/username/data"

            rapid_manager = RAPID()
            rapid_manager.generate_usgs_avg_daily_flows_opt(
                reach_id_gage_id_file=join(main_path, "usgsgage_id_comid.csv"),
                start_datetime=datetime.datetime(2000,1,1),
                end_datetime=datetime.datetime(2014,12,31),
                out_streamflow_file=join(main_path,"streamflow_2000_2014.csv"),
                out_stream_id_file=join(main_path,"streamid_2000_2014.csv")
            )


        Complex Example:

        .. code:: python

            import datetime
            from os.path import join
            from RAPIDpy import RAPID

            main_path = "/home/username/data"

            rapid_manager = RAPID(
                rapid_executable_location='~/work/rapid/run/rapid'
                use_all_processors=True,
                ZS_TauR=24*3600,
                ZS_dtR=15*60,
                ZS_TauM=365*24*3600,
                ZS_dtM=24*3600
            )

            rapid_manager.update_parameters(
                rapid_connect_file='../rapid-io/input/rapid_connect.csv',
                Vlat_file='../rapid-io/input/m3_riv.nc',
                riv_bas_id_file='../rapid-io/input/riv_bas_id.csv',
                k_file='../rapid-io/input/k.csv',
                x_file='../rapid-io/input/x.csv',
                Qout_file='../rapid-io/output/Qout.nc',
            )

            rapid_manager.update_reach_number_data()
            rapid_manager.update_simulation_runtime()
            rapid_manager.generate_usgs_avg_daily_flows_opt(
                reach_id_gage_id_file=join(main_path, "usgsgage_id_comid.csv"),
                start_datetime=datetime.datetime(2000,1,1),
                end_datetime=datetime.datetime(2014,12,31),
                out_streamflow_file=join(main_path,"streamflow_2000_2014.csv"),
                out_stream_id_file=join(main_path,"streamid_2000_2014.csv")
                )
            rapid_manager.run()

        """
        log("Generating avg streamflow file and stream id file "
            "required for calibration ...",
            "INFO")
        log("Generating avg streamflow file and stream id file "
            "required for calibration ...",
            "INFO")
        reach_id_gage_id_list = csv_to_list(reach_id_gage_id_file)
        gage_data_matrix = []
        valid_comid_list = []

        # add extra day as it includes the start date
        # (e.g. 7-5 is 2 days, but have data for 5,6,7, so +1)
        num_days_needed = (end_datetime-start_datetime).days + 1

        gage_id_list = []
        for row in reach_id_gage_id_list[1:]:
            station_id = row[1]
            if len(row[1]) == 7:
                station_id = '0' + row[1]
            gage_id_list.append(station_id)

        num_gage_id_list = np.array(gage_id_list, dtype=np.int32)
        log("Querying Server for Data ...",
            "INFO")

        query_params = {
                        'format': 'json',
                        'sites': ",".join(gage_id_list),
                        'startDT': start_datetime.strftime("%Y-%m-%d"),
                        'endDT': end_datetime.strftime("%Y-%m-%d"),
                        'parameterCd': '00060',  # streamflow
                        'statCd': '00003'  # average
                       }
        response = get("http://waterservices.usgs.gov/nwis/dv",
                       params=query_params)

        if not response.ok:
            log("USGS query error ...",
                "WARNING")
            return

        requested_data = None
        try:
            requested_data = response.json()['value']['timeSeries']
        except IndexError:
            pass

        if requested_data is not None:
            for time_series in enumerate(requested_data):
                usgs_station_full_name = time_series[1]['name']
                usgs_station_id = usgs_station_full_name.split(":")[1]
                gage_data = []
                for time_step in time_series[1]['values'][0]['value']:
                    local_datetime = parse(time_step['dateTime'])
                    if local_datetime > end_datetime:
                        break

                    if local_datetime >= start_datetime:
                        if not time_step['value']:
                            log("MISSING DATA for USGS Station {0} {1} {2}"
                                .format(usgs_station_id,
                                        local_datetime,
                                        time_step['value']),
                                "WARNING")
                        gage_data.append(
                            float(time_step['value']) / 35.3146667)

                try:
                    # get where streamids associated with USGS station ID
                    streamid_index = \
                        np.where(num_gage_id_list ==
                                 int(float(usgs_station_id)))[0][0]+1
                except (IndexError, ValueError):
                    log("USGS Station {0} not found in list ..."
                        .format(usgs_station_id),
                        "WARNING")
                    raise

                if len(gage_data) == num_days_needed:
                    gage_data_matrix.append(gage_data)
                    valid_comid_list.append(
                        reach_id_gage_id_list[streamid_index][0])
                else:
                    log("StreamID {0} USGS Station {1} MISSING {2} "
                        "DATA VALUES".format(
                            reach_id_gage_id_list[streamid_index][0],
                            usgs_station_id,
                            num_days_needed-len(gage_data)),
                        "WARNING")

        if gage_data_matrix and valid_comid_list:
            log("Writing Output ...",
                "INFO")
            np_array = np.array(gage_data_matrix).transpose()
            with open_csv(out_streamflow_file, 'w') as gage_data:
                wgd = csvwriter(gage_data)
                for row in np_array:
                    wgd.writerow(row)

            with open_csv(out_stream_id_file, 'w') as comid_data:
                wcd = csvwriter(comid_data)
                for row in valid_comid_list:
                    wcd.writerow([int(float(row))])

            # set parameters for RAPID run
            self.IS_obs_tot = len(valid_comid_list)
            self.obs_tot_id_file = out_stream_id_file
            self.Qobs_file = out_streamflow_file
            self.IS_obs_use = len(valid_comid_list)
            self.obs_use_id_file = out_stream_id_file
        else:
            log("No valid data returned ...",
                "WARNING")

'''-------------------------------------------------------------------------------
 Source Name: CreateInflowFileFromWRFHydroRunoff.py
 Author:      Environmental Systems Research Institute Inc.
 Updated by:  Alan D. Snow, US Army ERDC
 Description: Creates RAPID inflow file based on the WRF_Hydro land model output
              and the weight table previously created.
 History:     Initial coding - 10/17/2014, version 1.0
 ------------------------------------------------------------------------------'''
from .CreateInflowFileFromLDASRunoff import CreateInflowFileFromLDASRunoff

class CreateInflowFileFromWRFHydroRunoff(CreateInflowFileFromLDASRunoff):
    def __init__(self, lat_dim="south_north",
                 lon_dim="west_east",
                 lat_var="XLAT",
                 lon_var="XLONG",
                 surface_runoff_var="SFROFF",
                 subsurface_runoff_var="UDROFF",
                 time_step_seconds=1*3600):
                     
        """Define the tool (tool name is the name of the class)."""
        
        super(CreateInflowFileFromWRFHydroRunoff, self).__init__(lat_dim, lon_dim, lat_var, lon_var, 
                                                                 surface_runoff_var, subsurface_runoff_var, 
                                                                 time_step_seconds)
                                                                 
        self.header_wt = ['rivid', 'area_sqm', 'west_east', 'south_north', 'npoints']
        # According to David Gochis, underground runoff is "a major fraction of total river flow in most places"
        self.dims_oi = ['Time', lat_dim, lon_dim]

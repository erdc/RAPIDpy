"""
   CreateInflowFileFromLDASRunoff.py
   RAPIDpy

   Created by Alan D. Snow, 2016
   Adapted from CreateInflowFileFromWRFHydroRunoff.py.
   License: BSD-3-Clause
"""
from .CreateInflowFileFromLDASRunoff import CreateInflowFileFromLDASRunoff


class CreateInflowFileFromWRFHydroRunoff(CreateInflowFileFromLDASRunoff):
    """Create Inflow File From WRF-Hydro Runoff

    Base class for creating RAPID NetCDF input
    of water inflow based on WRF-Hydro
    runoff and previously created weight table.

    According to David Gochis, underground runoff is
    "a major fraction of total river flow in most places"
    """
    land_surface_model_name = "WRF-Hydro"
    header_wt = ['rivid', 'area_sqm', 'west_east', 'south_north', 'npoints']

    def __init__(self, lat_dim="south_north",
                 lon_dim="west_east",
                 lat_var="XLAT",
                 lon_var="XLONG",
                 surface_runoff_var="SFROFF",
                 subsurface_runoff_var="UDROFF"):

        """Define the tool (tool name is the name of the class)."""
        self.dims_oi = ['Time', lat_dim, lon_dim]

        super(CreateInflowFileFromWRFHydroRunoff, self).\
            __init__(lat_dim, lon_dim, lat_var, lon_var,
                     [surface_runoff_var, subsurface_runoff_var])

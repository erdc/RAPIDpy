GIS Tools
=========

For these tools to work, you need GIS dependencies installed (See :doc:`installation`).

There are also tools by Esri for ArcMap located here:

- https://github.com/Esri/python-toolbox-for-rapid
- https://github.com/erdc-cm/python-toolbox-for-rapid

To generate your own network from a DEM see :doc:`gis_stream_network`

Static RAPID Files
------------------

To prepare the static RAPID files (rapid\_connect.csv, riv\_bas\_id.csv,
kfac.csv, k.csv, x.csv, comid\_lat\_lon\_z.csv)

.. code:: python

    from RAPIDpy.gis.workflow import CreateAllStaticRAPIDFiles
    if __name__=="__main__":
        in_drainage_line = "/path/to/drainage_line.shp"
        in_catchment = "/path/to/catchment.shp"
        rapid_folder = "/path/to/rapid/output"
        CreateAllStaticRAPIDFiles(in_drainage_line=in_drainage_line,
                                  river_id="HydroID",
                                  length_id="LENGTHKM",
                                  slope_id="SLOPE",
                                  next_down_river_id="NextDownID",
                                  rapid_output_folder=rapid_folder,
                                  kfac_celerity=1000.0/3600.0, #this is the default (1km/hr converted to m/s)
                                  kfac_formula_type=3, #default is 3
                                  kfac_length_units="km", #length units (default is km)
                                  lambda_k=0.35, #k multiplying factor (default is 0.35)
                                  x_value=0.3, #Muskingum X [0-0.5] (default is 0.3)
                                  nhdplus=False, #set to true if you have VAA atributes added (FROMNODE, TONODE, DIVERGENCE)
                                 )

Weight Table Files
------------------

The weight tables are generated using an area weighted method based on
Esri's RAPID\_Toolbox.

To generate the ECMWF weight table files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    from RAPIDpy.gis.workflow import CreateAllStaticECMWFFiles
    if __name__=="__main__":
        in_catchment = "/path/to/catchment.shp"
        rapid_folder = "/path/to/rapid/output"
        rapid_connect = "/path/to/rapid_connect.csv"
        CreateAllStaticECMWFFiles(in_catchment=in_catchment,
                                  catchment_river_id="HydroID",
                                  rapid_output_folder=rapid_folder,
                                  rapid_connect_file=rapid_connect,
                                  )

To generate the NLDAS/GLDAS weight table files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    from RAPIDpy.gis.workflow import CreateWeightTableLDAS
    if __name__=="__main__":
        in_ldas_nc = "/path/to/nldas.nc"
        in_catchment = "/path/to/catchment.shp"
        rapid_connect = "/path/to/rapid_connect.csv"
        weight_file = "/path/to/rapid/output/weight_nldas.csv"
        CreateWeightTableLDAS(in_ldas_nc=in_ldas_nc,
                              in_nc_lon_var="lon_110",
                              in_nc_lat_var="lat_110",
                              in_catchment_shapefile=in_catchment, 
                              river_id="HydroID",
                              in_rapid_connect=rapid_connect, 
                              out_weight_table=weight_file)

How it works:
-------------

Snow, Alan D., Scott D. Christensen, Nathan R. Swain, E. James Nelson,
Daniel P. Ames, Norman L. Jones, Deng Ding, Nawajish S. Noman, Cedric H.
David, Florian Pappenberger, and Ervin Zsoter, 2016. A High-Resolution
National-Scale Hydrologic Forecast System from a Global Ensemble Land
Surface Model. *Journal of the American Water Resources Association
(JAWRA)* 1-15, DOI: 10.1111/1752-1688.12434

Snow, Alan Dee, "A New Global Forecasting Model to Produce
High-Resolution Stream Forecasts" (2015). All Theses and Dissertations.
Paper 5272. http://scholarsarchive.byu.edu/etd/5272

Calibration Tools
=================

Getting USGS Daily Gage Data
----------------------------

You can use this to download USGS daily avg flow data. BE CAREFUL
because overuse will get you blocked from downloading data.

Example reach\_id\_gage\_id\_file:

::

    COMID, USGS_GAGE_ID
    2000, 503944
    ...

Example:

.. code:: python

    from os.path import join
    main = "/home/username/data"
    rapid_manager.generate_usgs_avg_daily_flows_opt(reach_id_gage_id_file=join(main,"mississippi_usgsgage_id_comid.csv"),
                                                    start_datetime=datetime.datetime(2000,1,1),
                                                    end_datetime=datetime.datetime(2014,12,31),
                                                    out_streamflow_file=join(main,"streamflow_2000_2014.csv"), 
                                                    out_stream_id_file=join(main,"streamid_2000_2014.csv"))

Muskingum K File
----------------

Generate Muskingum Kfac file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To use this tool, you will need the GIS functions enabled (See:
https://github.com/erdc-cm/RAPIDpy/wiki/Installation).

There are three formula options:

1. River Length/Celerity
2. Eta\*River Length/Sqrt(River Slope)
3. Eta\*River Length/Sqrt(River Slope) - filters outliers using [0.05,
   0.95]

.. code:: python

    from RAPIDpy.gis.muskingum import CreateMuskingumKfacFile
    if __name__=="__main__":
        in_drainage_line = "/path/to/drainage_line.shp"
        rapid_connect_file = "/path/to/rapid/input/rapid_connect.csv"
        out_kfac_file = "/path/to/rapid/input/kfac.csv"
        CreateMuskingumKfacFile(in_drainage_line,
                                stream_id="HydroID",
                                length_id="LENGTHKM",
                                slope_id="SLOPE",
                                celerity=1000.0/3600.0, #units are m/s
                                formula_type=3, #default is 3
                                in_connectivity_file=rapid_connect_file,
                                out_kfac_file=out_kfac_file,
                                length_units="km", #length units (default is km))

Generate Muskingum K file
~~~~~~~~~~~~~~~~~~~~~~~~~

After calibrating with RAPID and you get a value for lambda, use this to
generate the Muskigum K file.

.. code:: python

    from RAPIDpy.gis.muskingum import CreateMuskingumKfacFile
    if __name__=="__main__":
        in_kfac_file = "/path/to/rapid/input/kfac.csv"
        out_k_file = "/path/to/rapid/input/k.csv"
        CreateMuskingumKFile(lambda_k=0.35,
                             in_kfac_file=in_kfac_file,
                             out_k_file=out_k_file)

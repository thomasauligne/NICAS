#!/bin/ksh
# ----------------------------------------------------------------------
# Korn shell script: gfs/links.ksh
# Author: Benjamin Menetrier
# Licensing: this code is distributed under the CeCILL-C license
# Copyright © 2017 METEO-FRANCE
# ----------------------------------------------------------------------

# Generate grid with ncks
ORIGIN_FILE=../../../../data/GFS/sfg_2014040100_fhr06s_mem001.nc4
rm -f grid.nc
ncks -O -v latitude,longitude,level,ak,bk ${ORIGIN_FILE} grid.nc

# Generate grid_SCRIP.nc with NCL
rm -f grid_SCRIP.nc
cat<<EOFNAM >ncl_request.ncl
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

data = addfile("grid.nc","r")
longitude = data->longitude
latitude = data->latitude

opt = True
opt@PrintTimings = True
opt@ForceOverwrite = True
rectilinear_to_SCRIP("grid_SCRIP.nc",latitude,longitude,opt)

end
EOFNAM
ncl ncl_request.ncl
rm -f ncl_request.ncl

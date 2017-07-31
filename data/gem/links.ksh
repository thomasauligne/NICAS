#!/bin/ksh
# ----------------------------------------------------------------------
# Korn shell script: gem/links.ksh
# Author: Benjamin Menetrier
# Licensing: this code is distributed under the CeCILL-C license
# Copyright © 2017 METEO-FRANCE
# ----------------------------------------------------------------------

# Generate grid with ncks
ORIGIN_FILE=../../../../data/GEM/2014101706_006_0001.nc
rm -f grid.nc
ncks -O -v lat,lon,lev,ap,b ${ORIGIN_FILE} grid.nc

# Generate grid_SCRIP.nc with NCL
rm -f grid_SCRIP.nc
cat<<EOFNAM >ncl_request.ncl
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

data = addfile("grid.nc","r")
lon = data->lon
lat = data->lat

opt = True
opt@PrintTimings = True
opt@ForceOverwrite = True
rectilinear_to_SCRIP("grid_SCRIP.nc",lat,lon,opt)

end
EOFNAM
ncl ncl_request.ncl
rm -f ncl_request.ncl

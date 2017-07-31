#!/bin/ksh
# ----------------------------------------------------------------------
# Korn shell script: nemo/links.ksh
# Author: Benjamin Menetrier
# Licensing: this code is distributed under the CeCILL-C license
# Copyright © 2017 METEO-FRANCE
# ----------------------------------------------------------------------

# Generate grid.nc with ncks
ORIGIN_FILE=../../../../data/NEMO/ENSEMBLES/mesh_mask
rm -f grid.nc
ncks -O -v nav_lat,nav_lon,tmask,e1t,e2t ${ORIGIN_FILE} grid.nc

# Generate grid_SCRIP.nc with NCL
rm -f grid_SCRIP.nc
cat<<EOFNAM >ncl_request.ncl
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

data = addfile("grid.nc","r")
nav_lon = data->nav_lon
nav_lat = data->nav_lat

opt = True
opt@PrintTimings = True
opt@ForceOverwrite = True
curvilinear_to_SCRIP("grid_SCRIP.nc",nav_lat,nav_lon,opt)

end
EOFNAM
ncl ncl_request.ncl
rm -f ncl_request.ncl

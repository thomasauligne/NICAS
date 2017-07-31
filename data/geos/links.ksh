#!/bin/ksh
# ----------------------------------------------------------------------
# Korn shell script: geos/links.ksh
# Author: Benjamin Menetrier
# Licensing: this code is distributed under the CeCILL-C license
# Copyright © 2017 METEO-FRANCE
# ----------------------------------------------------------------------

# Generate grid with ncks and ncwa
ORIGIN_FILE=../../../../data/GEOS/GEOS.fp.fcst.inst3_3d_asm_Nv.20160724_00+20160727_0000.V01.nc4
rm -f grid.nc
ncks -O -v lat,lon ${ORIGIN_FILE} grid.nc
ncwa -O -v PL -a time,lat,lon ${ORIGIN_FILE} pressure.nc
ncks -A -v PL pressure.nc grid.nc
rm -f pressure.nc

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

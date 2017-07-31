#!/bin/ksh
# ----------------------------------------------------------------------
# Korn shell script: wrf/links.ksh
# Author: Benjamin Menetrier
# Licensing: this code is distributed under the CeCILL-C license
# Copyright © 2017 METEO-FRANCE
# ----------------------------------------------------------------------

# Generate grid.nc with ncks and ncwa
ORIGIN_FILE=../../../../data/WRF/wrfda/2008020100/wrfout_d01_2008-02-01_00:00:00.nc
rm -f grid.nc
ncks -O -v XLONG,XLAT ${ORIGIN_FILE} grid.nc
ncwa -O -v PB -a Time,south_north,west_east ${ORIGIN_FILE} pressure.nc
ncks -A -v PB pressure.nc grid.nc
rm -f pressure.nc

# Generate grid_SCRIP.nc with NCL
rm -f grid_SCRIP.nc
cat<<EOFNAM >ncl_request.ncl
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

data = addfile("grid.nc","r")
XLONG = data->XLONG(0,:,:)
XLAT = data->XLAT(0,:,:)

opt = True
opt@PrintTimings = True
opt@ForceOverwrite = True
curvilinear_to_SCRIP("grid_SCRIP.nc",XLAT,XLONG,opt)

end
EOFNAM
ncl ncl_request.ncl
rm -f ncl_request.ncl

!----------------------------------------------------------------------
! Module: module_wrf.f90
!> Purpose: WRF model routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module model_wrf

use module_namelist, only: datadir,levs,logpres
use netcdf
use tools_const, only: deg2rad,req
use tools_missing, only: msvalr,msr,isanynotmsr
use tools_nc, only: ncerr,ncfloat
use type_sdata, only: sdatatype,sdata_alloc
implicit none

private
public :: model_wrf_coord,model_wrf_read,model_wrf_write

contains

!----------------------------------------------------------------------
! Subroutine: model_wrf_coord
!> Purpose: get WRF coordinates
!----------------------------------------------------------------------
subroutine model_wrf_coord(sdata)

implicit none

! Passed variables
type(sdatatype),intent(inout) :: sdata !< Sampling data

! Local variables
integer :: ic0
integer :: ncid,nlon_id,nlat_id,nlev_id,lon_id,lat_id,pres_id
real(kind=4) :: dx,dy
real(kind=4),allocatable :: lon(:,:),lat(:,:),pres(:)
character(len=1024) :: subr = 'model_wrf_coord'

! Open file and get dimensions
call ncerr(subr,nf90_open(trim(datadir)//'/grid.nc',nf90_nowrite,ncid))
call ncerr(subr,nf90_inq_dimid(ncid,'west_east',nlon_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nlon_id,len=sdata%nlon))
call ncerr(subr,nf90_inq_dimid(ncid,'south_north',nlat_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nlat_id,len=sdata%nlat))
sdata%nc0 = sdata%nlon*sdata%nlat
call ncerr(subr,nf90_inq_dimid(ncid,'bottom_top',nlev_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nlev_id,len=sdata%nlev))

! Allocation
allocate(lon(sdata%nlon,sdata%nlat))
allocate(lat(sdata%nlon,sdata%nlat))
allocate(pres(sdata%nlev))

! Read data and close file
call ncerr(subr,nf90_inq_varid(ncid,'XLONG',lon_id))
call ncerr(subr,nf90_inq_varid(ncid,'XLAT',lat_id))
call ncerr(subr,nf90_inq_varid(ncid,'PB',pres_id))
call ncerr(subr,nf90_get_var(ncid,lon_id,lon,(/1,1,1/),(/sdata%nlon,sdata%nlat,1/)))
call ncerr(subr,nf90_get_var(ncid,lat_id,lat,(/1,1,1/),(/sdata%nlon,sdata%nlat,1/)))
call ncerr(subr,nf90_get_var(ncid,pres_id,pres))
call ncerr(subr,nf90_get_att(ncid,nf90_global,'DX',dx))
call ncerr(subr,nf90_get_att(ncid,nf90_global,'DY',dy))
call ncerr(subr,nf90_close(ncid))

! Compute normalized area
allocate(sdata%area(sdata%nl0))
sdata%area = float(sdata%nlon*sdata%nlat)*dx*dy/req**2

! Conversion array
allocate(sdata%ic0(sdata%nc0))
do ic0=1,sdata%nc0
   sdata%ic0(ic0) = ic0
end do

! Conversion to radian
lon = lon*real(deg2rad,kind=4)
lat = lat*real(deg2rad,kind=4)

! Pack
call sdata_alloc(sdata)
sdata%lon = pack(real(lon,kind(1.0)),mask=.true.)
sdata%lat = pack(real(lat,kind(1.0)),mask=.true.)
sdata%mask = .true.

! Vertical unit
if (logpres) then
   sdata%vunit = log(pres(levs(1:sdata%nl0)))
else
   sdata%vunit = float(levs(1:sdata%nl0))
end if

! Release memory
deallocate(lon)
deallocate(lat)
deallocate(pres)

end subroutine model_wrf_coord

!----------------------------------------------------------------------
! Subroutine: model_wrf_read
!> Purpose: read WRF field
!----------------------------------------------------------------------
subroutine model_wrf_read(ncid,varname,sdata,fld)

implicit none

! Passed variables
integer,intent(in) :: ncid                   !< NetCDF file ID
character(len=*),intent(in) :: varname       !< Variable name
type(sdatatype),intent(in) :: sdata          !< Sampling data
real,intent(out) :: fld(sdata%nc0,sdata%nl0) !< Read field

! Local variables
integer :: il0
integer :: fld_id
real(kind=4),allocatable :: fld_loc(:,:)
character(len=1024) :: subr = 'model_wrf_read'

! Initialize field
call msr(fld)

! Allocation
allocate(fld_loc(sdata%nlon,sdata%nlat))

! Get variable id
call ncerr(subr,nf90_inq_varid(ncid,trim(varname),fld_id))

! Read variable
do il0=1,sdata%nl0
   call ncerr(subr,nf90_get_var(ncid,fld_id,fld_loc,(/1,1,levs(il0),1/),(/sdata%nlon,sdata%nlat,1,1/)))
   fld(:,il0) = pack(real(fld_loc,kind(1.0)),mask=.true.)
end do

! Release memory
deallocate(fld_loc)

end subroutine model_wrf_read


!----------------------------------------------------------------------
! Subroutine: model_wrf_write
!> Purpose: write WRF field
!----------------------------------------------------------------------
subroutine model_wrf_write(ncid,varname,sdata,fld)

implicit none

! Passed variables
integer,intent(in) :: ncid                  !< NetCDF file ID
character(len=*),intent(in) :: varname      !< Variable name
type(sdatatype),intent(in) :: sdata         !< Sampling data
real,intent(in) :: fld(sdata%nc0,sdata%nl0) !< Written field

! Local variables
integer :: il0,ierr
integer :: nlon_id,nlat_id,nlev_id,nt_id,fld_id
real :: fld_loc(sdata%nlon,sdata%nlat)
logical :: mask_unpack(sdata%nlon,sdata%nlat)
character(len=1024) :: subr = 'model_wrf_write'

! Initialization
mask_unpack = .true.

! Get variable id
ierr = nf90_inq_varid(ncid,trim(varname),fld_id)

! Define dimensions and variable if necessary
if (ierr/=nf90_noerr) then
   call ncerr(subr,nf90_redef(ncid))
   ierr = nf90_inq_dimid(ncid,'west_east',nlon_id)
   if (ierr/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'west_east',sdata%nlon,nlon_id))
   ierr = nf90_inq_dimid(ncid,'south_north',nlat_id)
   if (ierr/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'south_north',sdata%nlat,nlat_id))
   ierr = nf90_inq_dimid(ncid,'bottom_top',nlev_id)
   if (ierr/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'bottom_top',sdata%nl0,nlev_id))
   ierr = nf90_inq_dimid(ncid,'Time',nt_id)
   if (ierr/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'Time',1,nt_id))
   call ncerr(subr,nf90_def_var(ncid,trim(varname),ncfloat,(/nlon_id,nlat_id,nlev_id,nt_id/),fld_id))
   call ncerr(subr,nf90_put_att(ncid,fld_id,'_FillValue',msvalr))
   call ncerr(subr,nf90_enddef(ncid))
end if

! Write data
do il0=1,sdata%nl0
   if (isanynotmsr(fld(:,il0))) then
      call msr(fld_loc)
      fld_loc = unpack(fld(:,il0),mask=mask_unpack,field=fld_loc)
      call ncerr(subr,nf90_put_var(ncid,fld_id,fld_loc,(/1,1,il0,1/),(/sdata%nlon,sdata%nlat,1,1/)))
   end if
end do

end subroutine model_wrf_write

end module model_wrf

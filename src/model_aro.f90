!----------------------------------------------------------------------
! Module: module_aro.f90
!> Purpose: AROME model routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module model_aro

use module_namelist, only: datadir,levs,logpres
use netcdf
use tools_const, only: deg2rad,rad2deg,req,ps
use tools_missing, only: msvalr,msr,isanynotmsr
use tools_nc, only: ncerr,ncfloat
use type_sdata, only: sdatatype,sdata_alloc
implicit none

private
public :: model_aro_coord,model_aro_read,model_aro_write

contains

!----------------------------------------------------------------------
! Subroutine: model_aro_coord
!> Purpose: load AROME coordinates
!----------------------------------------------------------------------
subroutine model_aro_coord(sdata)

implicit none

! Passed variables
type(sdatatype),intent(inout) :: sdata !< Sampling data

! Local variables
integer :: il0,ic0
integer :: ncid,nlon_id,nlat_id,nlev_id,pp_id,lon_id,lat_id,cmask_id,a_id,b_id
real :: dx,dy
real(kind=8),allocatable :: lon(:,:),lat(:,:),cmask(:,:),a(:),b(:)
character(len=1024) :: subr = 'model_aro_coord'

! Open file and get dimensions
call ncerr(subr,nf90_open(trim(datadir)//'/grid.nc',nf90_nowrite,ncid))
call ncerr(subr,nf90_inq_dimid(ncid,'X',nlon_id))
call ncerr(subr,nf90_inq_dimid(ncid,'Y',nlat_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nlon_id,len=sdata%nlon))
call ncerr(subr,nf90_inquire_dimension(ncid,nlat_id,len=sdata%nlat))
sdata%nc0 = sdata%nlon*sdata%nlat
call ncerr(subr,nf90_inq_dimid(ncid,'Z',nlev_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nlev_id,len=sdata%nlev))

! Allocation
allocate(lon(sdata%nlon,sdata%nlat))
allocate(lat(sdata%nlon,sdata%nlat))
allocate(cmask(sdata%nlon,sdata%nlat))
allocate(a(sdata%nlev+1))
allocate(b(sdata%nlev+1))

! Read data and close file
call ncerr(subr,nf90_inq_varid(ncid,'longitude',lon_id))
call ncerr(subr,nf90_inq_varid(ncid,'latitude',lat_id))
call ncerr(subr,nf90_inq_varid(ncid,'cmask',cmask_id))
call ncerr(subr,nf90_inq_varid(ncid,'hybrid_coef_A',a_id))
call ncerr(subr,nf90_inq_varid(ncid,'hybrid_coef_B',b_id))
call ncerr(subr,nf90_inq_varid(ncid,'Projection_parameters',pp_id))
call ncerr(subr,nf90_get_var(ncid,lon_id,lon))
call ncerr(subr,nf90_get_var(ncid,lat_id,lat))
call ncerr(subr,nf90_get_var(ncid,cmask_id,cmask))
call ncerr(subr,nf90_get_var(ncid,a_id,a))
call ncerr(subr,nf90_get_var(ncid,b_id,b))
call ncerr(subr,nf90_get_att(ncid,pp_id,'x_resolution',dx))
call ncerr(subr,nf90_get_att(ncid,pp_id,'y_resolution',dy))
call ncerr(subr,nf90_close(ncid))

! Compute normalized area
allocate(sdata%area(sdata%nl0))
sdata%area = float(sdata%nlon*sdata%nlat)*dx*dy/req**2

! Conversion array
allocate(sdata%ic0(sdata%nc0))
do ic0=1,sdata%nc0
   sdata%ic0(ic0) = ic0
end do

! Convert to radian
lon = lon*real(deg2rad,kind=8)
lat = lat*real(deg2rad,kind=8)

! Pack
call sdata_alloc(sdata)
sdata%lon = pack(real(lon,kind(1.0)),mask=.true.)
sdata%lat = pack(real(lat,kind(1.0)),mask=.true.)
sdata%mask = .true.

! Vertical unit
if (logpres) then
   do il0=1,sdata%nl0
      if (levs(il0)<=sdata%nlev) then
         sdata%vunit(il0) = log(0.5*(a(levs(il0))+a(levs(il0)+1))+0.5*(b(levs(il0))+b(levs(il0)+1))*ps)
      else
         sdata%vunit(il0) = log(ps)
      end if
   end do
else
   sdata%vunit = float(levs(1:sdata%nl0))
end if

! Release memory
deallocate(lon)
deallocate(lat)
deallocate(cmask)
deallocate(a)
deallocate(b)

end subroutine model_aro_coord

!----------------------------------------------------------------------
! Subroutine: model_aro_read
!> Purpose: read AROME field
!----------------------------------------------------------------------
subroutine model_aro_read(ncid,varname,sdata,fld)

implicit none

! Passed variables
integer,intent(in) :: ncid                   !< NetCDF file ID
character(len=*),intent(in) :: varname       !< Variable name
type(sdatatype),intent(in) :: sdata          !< Sampling data
real,intent(out) :: fld(sdata%nc0,sdata%nl0) !< Read field

! Local variables
integer :: il0
integer :: fld_id
real :: fld_loc(sdata%nlon,sdata%nlat,sdata%nl0)
character(len=1024) :: subr = 'model_aro_read'

! Initialize field
call msr(fld)

! Get variable id
call ncerr(subr,nf90_inq_varid(ncid,trim(varname),fld_id))

! Read data
do il0=1,sdata%nl0
   call ncerr(subr,nf90_get_var(ncid,fld_id,fld_loc(:,:,il0),(/1,1,levs(il0)/),(/sdata%nlon,sdata%nlat,1/)))
end do

! Pack data
do il0=1,sdata%nl0
   fld(:,il0) = pack(real(fld_loc(:,:,il0),kind(1.0)),mask=.true.)
end do

end subroutine model_aro_read

!----------------------------------------------------------------------
! Subroutine: model_aro_write
!> Purpose: write AROME field
!----------------------------------------------------------------------
subroutine model_aro_write(ncid,varname,sdata,fld)

implicit none

! Passed variables
integer,intent(in) :: ncid                  !< NetCDF file ID
character(len=*),intent(in) :: varname      !< Variable name
type(sdatatype),intent(in) :: sdata         !< Sampling data
real,intent(in) :: fld(sdata%nc0,sdata%nl0) !< Written field

! Local variables
integer :: il0,ierr
integer :: nlon_id,nlat_id,nlev_id,fld_id,lon_id,lat_id
real :: fld_loc(sdata%nlon,sdata%nlat)
logical :: mask_unpack(sdata%nlon,sdata%nlat)
character(len=1024) :: subr = 'model_aro_write'

! Initialization
mask_unpack = .true.

   ! Get variable id
   ierr = nf90_inq_varid(ncid,trim(varname),fld_id)

   ! Define dimensions and variable if necessary
   if (ierr/=nf90_noerr) then
      call ncerr(subr,nf90_redef(ncid))
      ierr = nf90_inq_dimid(ncid,'X',nlon_id)
      if (ierr/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'X',sdata%nlon,nlon_id))
      ierr = nf90_inq_dimid(ncid,'Y',nlat_id)
      if (ierr/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'Y',sdata%nlat,nlat_id))
      ierr = nf90_inq_dimid(ncid,'Z',nlev_id)
      if (ierr/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'Z',sdata%nl0,nlev_id))
      call ncerr(subr,nf90_def_var(ncid,trim(varname),ncfloat,(/nlon_id,nlat_id,nlev_id/),fld_id))
      call ncerr(subr,nf90_put_att(ncid,fld_id,'_FillValue',msvalr))
      call ncerr(subr,nf90_enddef(ncid))
   end if

   ! Write data
   do il0=1,sdata%nl0
      if (isanynotmsr(fld(:,il0))) then
         call msr(fld_loc)
         fld_loc = unpack(fld(:,il0),mask=mask_unpack,field=fld_loc)
         call ncerr(subr,nf90_put_var(ncid,fld_id,fld_loc,(/1,1,il0/),(/sdata%nlon,sdata%nlat,1/)))
      end if
   end do

! Write coordinates
ierr = nf90_inq_varid(ncid,'longitude',lon_id)
if (ierr/=nf90_noerr) then
   call ncerr(subr,nf90_redef(ncid))
   ierr = nf90_inq_dimid(ncid,'X',nlon_id)
   if (ierr/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'X',sdata%nlon,nlon_id))
   ierr = nf90_inq_dimid(ncid,'Y',nlat_id)
   if (ierr/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'Y',sdata%nlat,nlat_id))
   call ncerr(subr,nf90_def_var(ncid,'longitude',ncfloat,(/nlon_id,nlat_id/),lon_id))
   call ncerr(subr,nf90_def_var(ncid,'latitude',ncfloat,(/nlon_id,nlat_id/),lat_id))
   call ncerr(subr,nf90_enddef(ncid))
   fld_loc = unpack(sdata%lon*rad2deg,mask=mask_unpack,field=fld_loc)
   call ncerr(subr,nf90_put_var(ncid,lon_id,fld_loc))
   fld_loc = unpack(sdata%lat*rad2deg,mask=mask_unpack,field=fld_loc)
   call ncerr(subr,nf90_put_var(ncid,lat_id,fld_loc))
end if

end subroutine model_aro_write

end module model_aro

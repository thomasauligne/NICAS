!----------------------------------------------------------------------
! Module: model_interface.f90
!> Purpose: model routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module model_interface

use model_aro, only: model_aro_coord,model_aro_read,model_aro_write
use model_arp, only: model_arp_coord,model_arp_read,model_arp_write
use model_gem, only: model_gem_coord,model_gem_read,model_gem_write
use model_geos, only: model_geos_coord,model_geos_read,model_geos_write
use model_gfs, only: model_gfs_coord,model_gfs_read,model_gfs_write
use model_ifs, only: model_ifs_coord,model_ifs_read,model_ifs_write
use model_mpas, only: model_mpas_coord,model_mpas_read,model_mpas_write
use model_nemo, only: model_nemo_coord,model_nemo_read,model_nemo_write
use model_wrf, only: model_wrf_coord,model_wrf_read,model_wrf_write
use module_namelist, only: datadir,model,nl,levs,logpres
use netcdf
use tools_display, only: msgerror
use tools_missing, only: msvalr,msr
use tools_nc, only: ncfloat,ncerr
use type_mpl, only: mpl,mpl_bcast
use type_sdata, only: sdatatype
implicit none

private
public :: model_coord,model_read,model_write

contains

!----------------------------------------------------------------------
! Subroutine: model_coord
!> Purpose: get coordinates
!----------------------------------------------------------------------
subroutine model_coord(sdata)

implicit none

! Passed variables
type(sdatatype),intent(inout) :: sdata !< Sampling data

! Local variables
integer :: il0
logical :: same_mask

! TODO: change that one day
sdata%nl0 = nl

! Select model
if (trim(model)=='aro') then
   call model_aro_coord(sdata)
elseif (trim(model)=='arp') then
   call model_arp_coord(sdata)
elseif (trim(model)=='gem') then
   call model_gem_coord(sdata)
elseif (trim(model)=='geos') then
   call model_geos_coord(sdata)
elseif (trim(model)=='gfs') then
   call model_gfs_coord(sdata)
elseif (trim(model)=='ifs') then
   call model_ifs_coord(sdata)
elseif (trim(model)=='mpas') then
   call model_mpas_coord(sdata)
elseif (trim(model)=='nemo') then
   call model_nemo_coord(sdata)
elseif (trim(model)=='wrf') then
   call model_wrf_coord(sdata)
else
   call msgerror('wrong model')
end if

! Check if the mask is the same for all levels
same_mask = .true.
do il0=2,sdata%nl0
   same_mask = same_mask.and.(all((sdata%mask(:,il0).and.sdata%mask(:,1)).or.(.not.sdata%mask(:,il0).and..not.sdata%mask(:,1))))
end do

! Define number of independent levels
if (same_mask) then
   sdata%nl0i = 1
else
   sdata%nl0i = sdata%nl0
end if
write(mpl%unit,'(a7,a,i3)') '','Number of independent levels: ',sdata%nl0i

end subroutine model_coord

!----------------------------------------------------------------------
! Subroutine: model_read
!> Purpose: read model field
!----------------------------------------------------------------------
subroutine model_read(filename,varname,sdata,fld)

implicit none

! Passed variables
character(len=*),intent(in) :: filename      !< File name
character(len=*),intent(in) :: varname       !< Variable name
type(sdatatype),intent(in) :: sdata          !< Sampling data
real,intent(out) :: fld(sdata%nc0,sdata%nl0) !< Read field

! Local variables
integer :: ncid
character(len=1024) :: subr = 'model_read'

! Open file
call ncerr(subr,nf90_open(trim(datadir)//'/'//trim(filename),nf90_nowrite,ncid))

! Select model
if (trim(model)=='aro') then
   call model_aro_read(ncid,varname,sdata,fld)
elseif (trim(model)=='arp') then
   call model_arp_read(ncid,varname,sdata,fld)
elseif (trim(model)=='gem') then
   call model_gem_read(ncid,varname,sdata,fld)
elseif (trim(model)=='geos') then
   call model_geos_read(ncid,varname,sdata,fld)
elseif (trim(model)=='gfs') then
   call model_gfs_read(ncid,varname,sdata,fld)
elseif (trim(model)=='ifs') then
   call model_ifs_read(ncid,varname,sdata,fld)
elseif (trim(model)=='mpas') then
   call model_mpas_read(ncid,varname,sdata,fld)
elseif (trim(model)=='nemo') then
   call model_nemo_read(ncid,varname,sdata,fld)
elseif (trim(model)=='wrf') then
   call model_wrf_read(ncid,varname,sdata,fld)
else
   call msgerror('wrong model')
end if

! Close file
call ncerr(subr,nf90_close(ncid))

end subroutine model_read

!----------------------------------------------------------------------
! Subroutine: model_write
!> Purpose: write model field
!----------------------------------------------------------------------
subroutine model_write(filename,varname,sdata,fld)

implicit none

! Passed variables
character(len=*),intent(in) :: filename     !< File name
character(len=*),intent(in) :: varname      !< Variable name
type(sdatatype),intent(in) :: sdata         !< Sampling data
real,intent(in) :: fld(sdata%nc0,sdata%nl0) !< Written field

! Local variables
integer :: ierr
integer :: ncid
character(len=1024) :: subr = 'model_write_grid'

if (mpl%main) then
   ! Check if the file exists
   ierr = nf90_create(trim(datadir)//'/'//trim(filename),or(nf90_noclobber,nf90_64bit_offset),ncid)
   if (ierr/=nf90_noerr) then
      call ncerr(subr,nf90_open(trim(datadir)//'/'//trim(filename),nf90_write,ncid))
      call ncerr(subr,nf90_redef(ncid))
      call ncerr(subr,nf90_put_att(ncid,nf90_global,'_FillValue',msvalr))
   end if
   call ncerr(subr,nf90_enddef(ncid))

   ! Select model
   if (trim(model)=='aro') then
      call model_aro_write(ncid,varname,sdata,fld)
   elseif (trim(model)=='arp') then
      call model_arp_write(ncid,varname,sdata,fld)
   elseif (trim(model)=='gem') then
      call model_gem_write(ncid,varname,sdata,fld)
   elseif (trim(model)=='geos') then
      call model_geos_write(ncid,varname,sdata,fld)
   elseif (trim(model)=='gfs') then
      call model_gfs_write(ncid,varname,sdata,fld)
   elseif (trim(model)=='ifs') then
      call model_ifs_write(ncid,varname,sdata,fld)
   elseif (trim(model)=='mpas') then
      call model_mpas_write(ncid,varname,sdata,fld)
   elseif (trim(model)=='nemo') then
      call model_nemo_write(ncid,varname,sdata,fld)
   elseif (trim(model)=='wrf') then
      call model_wrf_write(ncid,varname,sdata,fld)
   else
      call msgerror('wrong model')
   end if

   ! Close file
   call ncerr(subr,nf90_close(ncid))
end if

end subroutine model_write

end module model_interface

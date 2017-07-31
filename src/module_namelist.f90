!----------------------------------------------------------------------
! Module: module_namelist
!> Purpose: namelist parameters management
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module module_namelist

use netcdf, only: nf90_put_att,nf90_global
use omp_lib, only: omp_get_num_procs
use tools_display, only: black,err,wng,ddis,msgerror,msgwarning
use tools_missing, only: msr
use tools_nc, only: ncerr
use type_mpl, only: mpl,mpl_bcast
implicit none

! Namelist parameters maximum sizes
integer,parameter :: nlmax = 200   !< Maximum number of levels
integer,parameter :: ndirmax = 100 !< Maximum number of diracs

! general_param
character(len=1024) :: datadir     !< Data directory
character(len=1024) :: prefix      !< Files prefix
character(len=1024) :: model       !< Model name ('aro','arp', 'gfs', 'ifs','mpas', 'nemo' or 'wrf')
logical :: colorlog                !< Add colors to the log (for display on terminal)
integer :: nl                      !< Number of levels
integer :: levs(nlmax)             !< Levels
logical :: new_param               !< Compute new parameters (if false, read file)
logical :: new_mpi                 !< Compute new mpi splitting (if false, read file)
logical :: check_adjoints          !< Test adjoints
logical :: check_pos_def           !< Test positive definiteness
logical :: check_mpi               !< Test single proc/multi-procs equivalence
logical :: check_dirac             !< Test NICAS application on diracs
integer :: ndir                    !< Number of diracs
real :: dirlon(ndirmax)            !< Diracs longitudes
real :: dirlat(ndirmax)            !< Diracs latitudes

! sampling_param
logical :: sam_default_seed        !< Default seed for random numbers
logical :: mask_check              !< Check that interpolations do not cross mask boundaries
integer :: ntry                    !< Number of tries to get the most separated point for the zero-separation sampling
integer :: nrep                    !< Number of replacement to improve homogeneity of the zero-separation sampling
logical :: logpres                 !< Use pressure logarithm as vertical coordinate (model level if .false.)

! nicas_param
logical :: lsqrt                   !< Square-root formulation
character(len=1024) :: Lbh_file    !< Horizontal length-scale file
real :: Lbh(nlmax)                 !< Horizontal length-scale
character(len=1024) :: Lbv_file    !< Vertical length-scale file
real :: Lbv(nlmax)                 !< Vertical length-scale
real :: resol                      !< Resolution
logical :: network                 !< Network-base convolution calculation (distance-based if false)
integer :: nproc                   !< Number of tasks
integer :: mpicom                  !< Number of communication steps

! Namelist blocks
namelist/general_param/datadir,prefix,colorlog,model,nl,levs,new_param,new_mpi,check_adjoints,check_pos_def,check_mpi,check_dirac,ndir,dirlon,dirlat
namelist/sampling_param/sam_default_seed,mask_check,ntry,nrep,logpres
namelist/nicas_param/lsqrt,Lbh_file,Lbh,Lbv_file,Lbv,resol,network,nproc,mpicom

interface namncwrite_param
  module procedure namncwrite_integer
  module procedure namncwrite_integer_array
  module procedure namncwrite_real
  module procedure namncwrite_real_array
  module procedure namncwrite_logical
  module procedure namncwrite_string
end interface

private
public :: datadir,prefix,colorlog,model,nl,levs,new_param,new_mpi,check_adjoints,check_pos_def,check_mpi,check_dirac,ndir,dirlon,dirlat
public :: sam_default_seed,mask_check,ntry,nrep,logpres
public :: lsqrt,Lbh_file,Lbh,Lbv_file,Lbv,resol,network,nproc,mpicom
public :: namread,namncwrite

contains

!----------------------------------------------------------------------
! Subroutine: namread
!> Purpose: read and check namelist parameters
!----------------------------------------------------------------------
subroutine namread

implicit none

! Local variables
integer :: il,idir

! general_param default
datadir = ''
prefix = ''
colorlog = .false.
model = ''
nl = -1
levs = -1
new_param = .false.
new_mpi = .false.
check_adjoints = .false.
check_pos_def = .false.
check_mpi = .false.
check_dirac = .false.
ndir = -1
dirlon = -999.0
dirlat = -999.0

! sampling_param default
sam_default_seed = .false.
mask_check = .false.
ntry = -1
nrep = -1
logpres = .false.

! nicas_param default
lsqrt = .false.
Lbh_file = ''
Lbh = -1.0
Lbv_file = ''
Lbv = -1.0
resol = -1.0
network = .false.
nproc = -1
mpicom = -1

if (mpl%main) then
   ! Read namelist
   read(*,nml=general_param)
   read(*,nml=sampling_param)
   read(*,nml=nicas_param)
end if

! Broadcast parameters
call mpl_bcast(datadir,mpl%ioproc)
call mpl_bcast(prefix,mpl%ioproc)
call mpl_bcast(colorlog,mpl%ioproc)
call mpl_bcast(model,mpl%ioproc)
call mpl_bcast(nl,mpl%ioproc)
call mpl_bcast(levs,mpl%ioproc)
call mpl_bcast(new_param,mpl%ioproc)
call mpl_bcast(new_mpi,mpl%ioproc)
call mpl_bcast(check_adjoints,mpl%ioproc)
call mpl_bcast(check_pos_def,mpl%ioproc)
call mpl_bcast(check_mpi,mpl%ioproc)
call mpl_bcast(check_dirac,mpl%ioproc)
call mpl_bcast(ndir,mpl%ioproc)
call mpl_bcast(dirlon,mpl%ioproc)
call mpl_bcast(dirlat,mpl%ioproc)
call mpl_bcast(sam_default_seed,mpl%ioproc)
call mpl_bcast(mask_check,mpl%ioproc)
call mpl_bcast(ntry,mpl%ioproc)
call mpl_bcast(nrep,mpl%ioproc)
call mpl_bcast(logpres,mpl%ioproc)
call mpl_bcast(lsqrt,mpl%ioproc)
call mpl_bcast(Lbh_file,mpl%ioproc)
call mpl_bcast(Lbh,mpl%ioproc)
call mpl_bcast(Lbv_file,mpl%ioproc)
call mpl_bcast(Lbv,mpl%ioproc)
call mpl_bcast(resol,mpl%ioproc)
call mpl_bcast(network,mpl%ioproc)
call mpl_bcast(nproc,mpl%ioproc)
call mpl_bcast(mpicom,mpl%ioproc)

! Setup display colors
if (colorlog) then
   black = char(27)//'[0;0m'
   err = char(27)//'[0;37;41;1m'
   wng = char(27)//'[0;37;42;1m'
else
   black = ' '
   err = ' '
   wng = ' '
end if
ddis = 5

! Check general_param
if (trim(datadir)=='') call msgerror('datadir not specified')
if (trim(prefix)=='') call msgerror('prefix not specified')
select case (trim(model))
case ('aro','arp','gem','geos','gfs','ifs','mpas','nemo','wrf')
case default
   call msgerror('wrong model')
end select
if (nl<=0) call msgerror('nl should be positive')
do il=1,nl
   if (levs(il)<=0) call msgerror('levs should be positive')
   if (count(levs(1:nl)==levs(il))>1) call msgerror('redundant levels')
end do
if (new_param.and.(.not.new_mpi)) then
   call msgwarning('new parameters calculation implies new MPI splitting, resetting new_mpi to .true.')
   new_mpi = .true.
end if
if (check_dirac) then
   if (ndir<1) call msgerror('ndir should be positive')
   do idir=1,ndir
      if ((dirlon(idir)<-180.0).or.(dirlon(idir)>180.0)) call msgerror('dirac longitude should lie between -180 and 180')
      if ((dirlat(idir)<-90.0).or.(dirlat(idir)>90.0)) call msgerror('dirac latitude should lie between -90 and 90')
   end do
end if

! Check sampling_param
if (ntry<=0) call msgerror('ntry should be positive')
if (nrep<0) call msgerror('nrep should be non-negative')
if (logpres) then
   select case (trim(model))
   case ('aro','arp','gem','geos','gfs','mpas','wrf')
   case default
      call msgwarning('pressure logarithm vertical coordinate is not available for this model, resetting to model level index')
      logpres = .false.
   end select
end if

! Check nicas_param
if (trim(Lbh_file)=='') then
   Lbh(1:nl) = Lbh(levs(1:nl))
   do il=1,nl
      if (Lbh(il)<tiny(1.0)) call msgerror('Lbh should be positive')
   end do
end if
if (trim(Lbv_file)=='') then
   Lbv(1:nl) = Lbv(levs(1:nl))
   do il=1,nl
      if (Lbv(il)<tiny(1.0)) call msgerror('Lbv should be positive')
   end do
end if
if (resol<tiny(1.0)) call msgerror('resol should be positive')
if (nproc<0) call msgerror('nproc should be non-negative')
if ((mpicom/=1).and.(mpicom/=2)) call msgerror('mpicom should be 1 or 2')

end subroutine namread

!----------------------------------------------------------------------
! Subroutine: namncwrite
!> Purpose: write namelist parameters as NetCDF attributes
!----------------------------------------------------------------------
subroutine namncwrite(ncid)

implicit none

! Passed variables
integer,intent(in) :: ncid !< NetCDF file id

! general_param
call namncwrite_param(ncid,'general_param_datadir',trim(datadir))
call namncwrite_param(ncid,'general_param_prefix',trim(prefix))
call namncwrite_param(ncid,'general_param_colorlog',colorlog)
call namncwrite_param(ncid,'general_param_model',trim(model))
call namncwrite_param(ncid,'general_param_nl',nl)
call namncwrite_param(ncid,'general_param_levs',nl,levs)
call namncwrite_param(ncid,'general_param_new_param',new_param)
call namncwrite_param(ncid,'general_param_new_mpi',new_mpi)
call namncwrite_param(ncid,'general_param_check_adjoints',check_adjoints)
call namncwrite_param(ncid,'general_param_check_pos_def',check_pos_def)
call namncwrite_param(ncid,'general_param_check_mpi',check_mpi)
call namncwrite_param(ncid,'general_param_check_dirac',check_dirac)
call namncwrite_param(ncid,'general_param_ndir',ndir)
call namncwrite_param(ncid,'general_param_dirlon',ndir,dirlon)
call namncwrite_param(ncid,'general_param_dirlat',ndir,dirlat)

! sampling_param
call namncwrite_param(ncid,'sampling_param_sam_default_seed',sam_default_seed)
call namncwrite_param(ncid,'sampling_param_mask_check',mask_check)
call namncwrite_param(ncid,'sampling_param_ntry',ntry)
call namncwrite_param(ncid,'sampling_param_nrep',nrep)
call namncwrite_param(ncid,'sampling_param_logpres',logpres)

! nicas_param
call namncwrite_param(ncid,'nicas_param_lsqrt',lsqrt)
call namncwrite_param(ncid,'nicas_param_Lbh_file',trim(Lbh_file))
call namncwrite_param(ncid,'nicas_param_Lbh',nl,Lbh)
call namncwrite_param(ncid,'nicas_param_Lbv_file',trim(Lbv_file))
call namncwrite_param(ncid,'nicas_param_Lbv',nl,Lbv)
call namncwrite_param(ncid,'nicas_param_resol',resol)
call namncwrite_param(ncid,'nicas_param_network',network)
call namncwrite_param(ncid,'nicas_param_nproc',nproc)
call namncwrite_param(ncid,'nicas_param_mpicom',mpicom)

end subroutine namncwrite

!----------------------------------------------------------------------
! Subroutine: namncwrite_integer
!> Purpose: write namelist integer as NetCDF attribute
!----------------------------------------------------------------------
subroutine namncwrite_integer(ncid,varname,var)

implicit none

! Passed variables
integer,intent(in) :: ncid             !< NetCDF file id
character(len=*),intent(in) :: varname !< Variable name
integer,intent(in) :: var              !< Integer

! Local variables
character(len=1024) :: subr='namncwrite_integer'

! Write integer
call ncerr(subr,nf90_put_att(ncid,nf90_global,trim(varname),var))

end subroutine namncwrite_integer

!----------------------------------------------------------------------
! Subroutine: namncwrite_integer_array
!> Purpose: write namelist integer array as NetCDF attribute
!----------------------------------------------------------------------
subroutine namncwrite_integer_array(ncid,varname,n,var)

implicit none

! Passed variables
integer,intent(in) :: ncid             !< NetCDF file id
character(len=*),intent(in) :: varname !< Variable name
integer,intent(in) :: n                !< Integer array size
integer,intent(in) :: var(n)           !< Integer array

! Local variables
integer :: i
character(len=1024) :: str,fullstr
character(len=1024) :: subr='namncwrite_integer_array'

! Write integer array as a string
if (n>0) then
   write(fullstr,'(i3.3)') var(1)
   do i=2,n
      write(str,'(i3.3)') var(i)
      fullstr = trim(fullstr)//':'//trim(str)
   end do
   call ncerr(subr,nf90_put_att(ncid,nf90_global,trim(varname),trim(fullstr)))
end if

end subroutine namncwrite_integer_array

!----------------------------------------------------------------------
! Subroutine: namncwrite_real
!> Purpose: write namelist real as NetCDF attribute
!----------------------------------------------------------------------
subroutine namncwrite_real(ncid,varname,var)

implicit none

! Passed variables
integer,intent(in) :: ncid             !< NetCDF file id
character(len=*),intent(in) :: varname !< Variable name
real,intent(in) :: var                 !< Real

! Local variables
character(len=1024) :: subr='namncwrite_real'

! Write real
call ncerr(subr,nf90_put_att(ncid,nf90_global,trim(varname),var))

end subroutine namncwrite_real

!----------------------------------------------------------------------
! Subroutine: namncwrite_real_array
!> Purpose: write namelist real array as NetCDF attribute
!----------------------------------------------------------------------
subroutine namncwrite_real_array(ncid,varname,n,var)

implicit none

! Passed variables
integer,intent(in) :: ncid             !< NetCDF file id
character(len=*),intent(in) :: varname !< Variable name
integer,intent(in) :: n                !< Real array size
real,intent(in) :: var(n)              !< Real array

! Local variables
integer :: i
character(len=1024) :: str,fullstr
character(len=1024) :: subr='namncwrite_real_array'

! Write real array as a string
if (n>0) then
   write(fullstr,'(e10.3)') var(1)
   do i=2,n
      write(str,'(e10.3)') var(i)
      fullstr = trim(fullstr)//':'//trim(str)
   end do
   call ncerr(subr,nf90_put_att(ncid,nf90_global,trim(varname),trim(fullstr)))
end if

end subroutine namncwrite_real_array

!----------------------------------------------------------------------
! Subroutine: namncwrite_logical
!> Purpose: write namelist logical as NetCDF attribute
!----------------------------------------------------------------------
subroutine namncwrite_logical(ncid,varname,var)

implicit none

! Passed variables
integer,intent(in) :: ncid             !< NetCDF file id
character(len=*),intent(in) :: varname !< Variable name
logical,intent(in) :: var              !< Logical

! Local variables
character(len=1024) :: subr='namncwrite_logical'

! Write logical as a string
if (var) then
   call ncerr(subr,nf90_put_att(ncid,nf90_global,trim(varname),'.true.'))
else
   call ncerr(subr,nf90_put_att(ncid,nf90_global,trim(varname),'.false.'))
end if

end subroutine namncwrite_logical

!----------------------------------------------------------------------
! Subroutine: namncwrite_string
!> Purpose: write namelist string as NetCDF attribute
!----------------------------------------------------------------------
subroutine namncwrite_string(ncid,varname,var)

implicit none

! Passed variables
integer,intent(in) :: ncid             !< NetCDF file id
character(len=*),intent(in) :: varname !< Variable name
character(len=*),intent(in) :: var     !< String

! Local variables
character(len=1024) :: subr='namncwrite_string'

! Write string
call ncerr(subr,nf90_put_att(ncid,nf90_global,trim(varname),trim(var)))

end subroutine namncwrite_string

end module module_namelist

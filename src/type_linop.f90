!----------------------------------------------------------------------
! Module: type_linop
!> Purpose: linear operator derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module type_linop

use netcdf
use omp_lib
use tools_display, only: msgerror
use tools_missing, only: msi,msr,isnotmsr
use tools_nc, only: ncfloat,ncerr
use type_mpl, only: mpl
implicit none

! Linear operator derived type
type linoptype
   character(len=1024) :: prefix !< Operator prefix (for I/O)
   integer :: n_src              !< Source vector size
   integer :: n_dst              !< Destination vector size
   integer :: n_s                !< Operator size
   integer,allocatable :: row(:) !< Output indices
   integer,allocatable :: col(:) !< Input indices
   real,allocatable :: S(:)      !< Coefficients
end type linoptype

interface linop_read
  module procedure linop_read_single
  module procedure linop_read_array
end interface
interface linop_write
  module procedure linop_write_single
  module procedure linop_write_array
end interface

logical :: check_data = .false. !< Activate data check for all linear operations

private
public :: linoptype
public :: linop_alloc,linop_dealloc,apply_linop,apply_linop_ad,apply_linop_sym,linop_read,linop_write

contains

!----------------------------------------------------------------------
! Subroutine: linop_alloc
!> Purpose: linear operator object allocation
!----------------------------------------------------------------------
subroutine linop_alloc(linop)

implicit none

! Passed variables
type(linoptype),intent(inout) :: linop !< Linear operator

! Allocation
allocate(linop%row(linop%n_s))
allocate(linop%col(linop%n_s))
allocate(linop%S(linop%n_s))

! Initialization
call msi(linop%row)
call msi(linop%col)
call msr(linop%S)

end subroutine linop_alloc

!----------------------------------------------------------------------
! Subroutine: linop_dealloc
!> Purpose: linear operator object deallocation
!----------------------------------------------------------------------
subroutine linop_dealloc(linop)

implicit none

! Passed variables
type(linoptype),intent(inout) :: linop !< Linear operator

! Release memory
deallocate(linop%row)
deallocate(linop%col)
deallocate(linop%S)

end subroutine linop_dealloc

!----------------------------------------------------------------------
! Subroutine: apply_linop
!> Purpose: apply linear operator
!----------------------------------------------------------------------
subroutine apply_linop(linop,fld_src,fld_dst)

implicit none

! Passed variables
type(linoptype),intent(in) :: linop      !< Linear operator
real,intent(in) :: fld_src(linop%n_src)  !< Source vector
real,intent(out) :: fld_dst(linop%n_dst) !< Destination vector

! Local variables
integer :: i_s

if (check_data) then
   ! Check linear operation
   if (minval(linop%col)<1) call msgerror('col<1 for linear operation '//trim(linop%prefix))
   if (maxval(linop%col)>linop%n_src) call msgerror('col>n_src for linear operation '//trim(linop%prefix))
   if (minval(linop%row)<1) call msgerror('row<1 for linear operation '//trim(linop%prefix))
   if (maxval(linop%row)>linop%n_dst) call msgerror('row>n_dst for linear operation '//trim(linop%prefix))
   if (any(isnan(linop%S))) call msgerror('NaN in S for linear operation '//trim(linop%prefix))

   ! Check input
   if (any(fld_src>huge(1.0))) call msgerror('Overflowing number in fld_src for linear operation '//trim(linop%prefix))
   if (any(isnan(fld_src))) call msgerror('NaN in fld_src for linear operation '//trim(linop%prefix))
end if

! Initialization
fld_dst = 0.0

! Apply weights
do i_s=1,linop%n_s
   fld_dst(linop%row(i_s)) = fld_dst(linop%row(i_s))+linop%S(i_s)*fld_src(linop%col(i_s))
end do

if (check_data) then
   ! Check output
   if (any(isnan(fld_dst))) call msgerror('NaN in fld_dst for linear operation '//trim(linop%prefix))
end if

end subroutine apply_linop

!----------------------------------------------------------------------
! Subroutine: apply_linop_ad
!> Purpose: apply linear operator, adjoint
!----------------------------------------------------------------------
subroutine apply_linop_ad(linop,fld_dst,fld_src)

implicit none

! Passed variables
type(linoptype),intent(in) :: linop      !< Linear operator
real,intent(in) :: fld_dst(linop%n_dst)  !< Destination vector
real,intent(out) :: fld_src(linop%n_src) !< Source vector

! Local variables
integer :: i_s

if (check_data) then
   ! Check linear operation
   if (minval(linop%col)<1) call msgerror('col<1 for adjoint linear operation '//trim(linop%prefix))
   if (maxval(linop%col)>linop%n_src) call msgerror('col>n_src for adjoint linear operation '//trim(linop%prefix))
   if (minval(linop%row)<1) call msgerror('row<1 for adjoint linear operation '//trim(linop%prefix))
   if (maxval(linop%row)>linop%n_dst) call msgerror('row>n_dst for adjoint linear operation '//trim(linop%prefix))
   if (any(isnan(linop%S))) call msgerror('NaN in S for adjoint linear operation '//trim(linop%prefix))

   ! Check input
   if (any(fld_dst>huge(1.0))) call msgerror('Overflowing number in fld_dst for adjoint linear operation '//trim(linop%prefix))
   if (any(isnan(fld_dst))) call msgerror('NaN in fld_dst for adjoint linear operation '//trim(linop%prefix))
end if

! Initialization
fld_src = 0.0

! Apply weights
do i_s=1,linop%n_s
   fld_src(linop%col(i_s)) = fld_src(linop%col(i_s))+linop%S(i_s)*fld_dst(linop%row(i_s))
end do

if (check_data) then
   ! Check output
   if (any(isnan(fld_src))) call msgerror('NaN in fld_src for adjoint linear operation '//trim(linop%prefix))
end if

end subroutine apply_linop_ad

!----------------------------------------------------------------------
! Subroutine: apply_linop_sym
!> Purpose: apply linear operator, symmetric
!----------------------------------------------------------------------
subroutine apply_linop_sym(linop,fld)

implicit none

! Passed variables
type(linoptype),intent(in) :: linop      !< Linear operator
real,intent(inout) :: fld(linop%n_src)   !< Source/destination vector

! Local variables
integer :: i_s,ithread
real :: fld_tmp(linop%n_dst,mpl%nthread)

if (check_data) then
   ! Check linear operation
   if (minval(linop%col)<1) call msgerror('col<1 for symmetric linear operation '//trim(linop%prefix))
   if (maxval(linop%col)>linop%n_src) call msgerror('col>n_src for symmetric linear operation '//trim(linop%prefix))
   if (minval(linop%row)<1) call msgerror('row<1 for symmetric linear operation '//trim(linop%prefix))
   if (maxval(linop%row)>linop%n_src) call msgerror('row>n_dst for symmetric linear operation '//trim(linop%prefix))
   if (any(isnan(linop%S))) call msgerror('NaN in S for symmetric linear operation '//trim(linop%prefix))

   ! Check input
   if (any(fld>huge(1.0))) call msgerror('Overflowing number in fld for symmetric linear operation '//trim(linop%prefix))
   if (any(isnan(fld))) call msgerror('NaN in fld for symmetric linear operation '//trim(linop%prefix))
end if

! Initialization
fld_tmp = 0.0

! Apply weights
!$omp parallel do private(i_s,ithread)
do i_s=1,linop%n_s
   ithread = omp_get_thread_num()+1
   fld_tmp(linop%row(i_s),ithread) = fld_tmp(linop%row(i_s),ithread)+linop%S(i_s)*fld(linop%col(i_s))
   fld_tmp(linop%col(i_s),ithread) = fld_tmp(linop%col(i_s),ithread)+linop%S(i_s)*fld(linop%row(i_s))
end do
!$omp end parallel do

! Sum over threads
do ithread=1,mpl%nthread
   fld = fld+fld_tmp(:,ithread)
end do

if (check_data) then
   ! Check output
   if (any(isnan(fld))) call msgerror('NaN in fld for symmetric linear operation '//trim(linop%prefix))
end if

end subroutine apply_linop_sym

!----------------------------------------------------------------------
! Subroutine: linop_read_single
!> Purpose: read single linear operator from a NetCDF file
!----------------------------------------------------------------------
subroutine linop_read_single(ncid,prefix,linop)

implicit none

! Passed variables
integer,intent(in) :: ncid             !< NetCDF file id
character(len=*),intent(in) :: prefix  !< Linear operator prefix
type(linoptype),intent(inout) :: linop !< Linear operator

! Local variables
integer :: info
integer :: n_s_id,row_id,col_id,S_id
character(len=1024) :: subr = 'linop_read_single'

! Copy prefix
linop%prefix = prefix

! Get operator size
info = nf90_inq_dimid(ncid,trim(prefix)//'_n_s',n_s_id)
if (info==nf90_noerr) then
   call ncerr(subr,nf90_inquire_dimension(ncid,n_s_id,len=linop%n_s))
else
   linop%n_s = 0
end if

if (linop%n_s>0) then
   ! Allocation
   call linop_alloc(linop)

   ! Get variables id
   call ncerr(subr,nf90_inq_varid(ncid,trim(prefix)//'_row',row_id))
   call ncerr(subr,nf90_inq_varid(ncid,trim(prefix)//'_col',col_id))
   call ncerr(subr,nf90_inq_varid(ncid,trim(prefix)//'_S',S_id))

   ! Get source/destination dimensions
   call ncerr(subr,nf90_get_att(ncid,nf90_global,trim(prefix)//'_n_src',linop%n_src))
   call ncerr(subr,nf90_get_att(ncid,nf90_global,trim(prefix)//'_n_dst',linop%n_dst))

   ! Get variables
   call ncerr(subr,nf90_get_var(ncid,row_id,linop%row))
   call ncerr(subr,nf90_get_var(ncid,col_id,linop%col))
   call ncerr(subr,nf90_get_var(ncid,S_id,linop%S))
end if

end subroutine linop_read_single

!----------------------------------------------------------------------
! Subroutine: linop_read_array
!> Purpose: read array of linear operators from a NetCDF file
!----------------------------------------------------------------------
subroutine linop_read_array(ncid,prefix,linop)

implicit none

! Passed variables
integer,intent(in) :: ncid                            !< NetCDF file id
character(len=*),intent(in) :: prefix                 !< Linear operator prefix
type(linoptype),allocatable,intent(inout) :: linop(:) !< Linear operators

! Local variables
integer :: info,narr,n_s_max,iarr
integer :: n_s_max_id,narr_id,n_s_id,row_id,col_id,S_id
character(len=1024) :: subr = 'linop_read_array'

! Get maximum operator size
info = nf90_inq_dimid(ncid,trim(prefix)//'_n_s_max',n_s_max_id)
if (info==nf90_noerr) then
   call ncerr(subr,nf90_inquire_dimension(ncid,n_s_max_id,len=n_s_max))
else
   n_s_max = 0
end if

! Get array size
info = nf90_inq_dimid(ncid,trim(prefix)//'_narr',narr_id)
if (info==nf90_noerr) then
   call ncerr(subr,nf90_inquire_dimension(ncid,narr_id,len=narr))
else
   narr = 0
end if

if ((narr>0).and.(n_s_max>0)) then
   ! Allocation
   allocate(linop(narr))

   ! Get variables id
   call ncerr(subr,nf90_inq_varid(ncid,trim(prefix)//'_n_s',n_s_id))
   call ncerr(subr,nf90_inq_varid(ncid,trim(prefix)//'_row',row_id))
   call ncerr(subr,nf90_inq_varid(ncid,trim(prefix)//'_col',col_id))
   call ncerr(subr,nf90_inq_varid(ncid,trim(prefix)//'_S',S_id))

   do iarr=1,narr
      ! Copy prefix
      linop(iarr)%prefix = prefix

      ! Get operator size
      call ncerr(subr,nf90_get_var(ncid,n_s_id,linop(iarr)%n_s,(/iarr/)))

      if (linop(iarr)%n_s>0) then
         ! Allocation
         call linop_alloc(linop(iarr))

         ! Get source/destination dimensions
         call ncerr(subr,nf90_get_att(ncid,nf90_global,trim(prefix)//'_n_src',linop(iarr)%n_src))
         call ncerr(subr,nf90_get_att(ncid,nf90_global,trim(prefix)//'_n_dst',linop(iarr)%n_dst))

         ! Get variables
         call ncerr(subr,nf90_get_var(ncid,row_id,linop(iarr)%row,(/1,iarr/),(/linop(iarr)%n_s,1/)))
         call ncerr(subr,nf90_get_var(ncid,col_id,linop(iarr)%col,(/1,iarr/),(/linop(iarr)%n_s,1/)))
         call ncerr(subr,nf90_get_var(ncid,S_id,linop(iarr)%S,(/1,iarr/),(/linop(iarr)%n_s,1/)))
      end if
   end do
end if

end subroutine linop_read_array

!----------------------------------------------------------------------
! Subroutine: linop_write_single
!> Purpose: write single linear operator to a NetCDF file
!----------------------------------------------------------------------
subroutine linop_write_single(ncid,linop)

implicit none

! Passed variables
integer,intent(in) :: ncid          !< NetCDF file id
type(linoptype),intent(in) :: linop !< Linear operator

! Local variables
integer :: n_s_id,row_id,col_id,S_id
character(len=1024) :: subr = 'linop_write_single'

if (linop%n_s>0) then
   ! Start definition mode
   call ncerr(subr,nf90_redef(ncid))

   ! Define dimensions
   call ncerr(subr,nf90_def_dim(ncid,trim(linop%prefix)//'_n_s',linop%n_s,n_s_id))

   ! Define variables
   call ncerr(subr,nf90_def_var(ncid,trim(linop%prefix)//'_row',nf90_int,(/n_s_id/),row_id))
   call ncerr(subr,nf90_def_var(ncid,trim(linop%prefix)//'_col',nf90_int,(/n_s_id/),col_id))
   call ncerr(subr,nf90_def_var(ncid,trim(linop%prefix)//'_S',ncfloat,(/n_s_id/),S_id))

   ! Write source/destination dimensions
   call ncerr(subr,nf90_put_att(ncid,nf90_global,trim(linop%prefix)//'_n_src',linop%n_src))
   call ncerr(subr,nf90_put_att(ncid,nf90_global,trim(linop%prefix)//'_n_dst',linop%n_dst))

   ! End definition mode
   call ncerr(subr,nf90_enddef(ncid))

   ! Put variables
   call ncerr(subr,nf90_put_var(ncid,row_id,linop%row))
   call ncerr(subr,nf90_put_var(ncid,col_id,linop%col))
   call ncerr(subr,nf90_put_var(ncid,S_id,linop%S))
end if

end subroutine linop_write_single

!----------------------------------------------------------------------
! Subroutine: linop_write_array
!> Purpose: write array of linear operators to a NetCDF file
!----------------------------------------------------------------------
subroutine linop_write_array(ncid,linop)

implicit none

! Passed variables
integer,intent(in) :: ncid             !< NetCDF file id
type(linoptype),intent(in) :: linop(:) !< Linear operator

! Local variables
integer :: narr,iarr,n_s_max
integer :: n_s_max_id,narr_id,n_s_id,row_id,col_id,S_id
character(len=1024) :: subr = 'linop_write_array'

! Array size
narr = size(linop)

! Maximum operator size
n_s_max = 0
do iarr=1,narr
   n_s_max = max(n_s_max,linop(iarr)%n_s)
end do

if ((narr>0).and.(n_s_max>0)) then
   ! Start definition mode
   call ncerr(subr,nf90_redef(ncid))

   ! Define dimension
   call ncerr(subr,nf90_def_dim(ncid,trim(linop(1)%prefix)//'_n_s_max',n_s_max,n_s_max_id))
   call ncerr(subr,nf90_def_dim(ncid,trim(linop(1)%prefix)//'_narr',narr,narr_id))

   ! Write source/destination dimensions
   call ncerr(subr,nf90_put_att(ncid,nf90_global,trim(linop(1)%prefix)//'_n_src',linop(1)%n_src))
   call ncerr(subr,nf90_put_att(ncid,nf90_global,trim(linop(1)%prefix)//'_n_dst',linop(1)%n_dst))

   ! Define variables
   call ncerr(subr,nf90_def_var(ncid,trim(linop(1)%prefix)//'_n_s',nf90_int,(/narr_id/),n_s_id))
   call ncerr(subr,nf90_def_var(ncid,trim(linop(1)%prefix)//'_row',nf90_int,(/n_s_max_id,narr_id/),row_id))
   call ncerr(subr,nf90_def_var(ncid,trim(linop(1)%prefix)//'_col',nf90_int,(/n_s_max_id,narr_id/),col_id))
   call ncerr(subr,nf90_def_var(ncid,trim(linop(1)%prefix)//'_S',ncfloat,(/n_s_max_id,narr_id/),S_id))

   ! End definition mode
   call ncerr(subr,nf90_enddef(ncid))

   do iarr=1,narr
      ! Put variables
      call ncerr(subr,nf90_put_var(ncid,n_s_id,linop(iarr)%n_s,(/iarr/)))
      call ncerr(subr,nf90_put_var(ncid,row_id,linop(iarr)%row,(/1,iarr/),(/linop(iarr)%n_s,1/)))
      call ncerr(subr,nf90_put_var(ncid,col_id,linop(iarr)%col,(/1,iarr/),(/linop(iarr)%n_s,1/)))
      call ncerr(subr,nf90_put_var(ncid,S_id,linop(iarr)%S,(/1,iarr/),(/linop(iarr)%n_s,1/)))
   end do
end if

end subroutine linop_write_array

end module type_linop

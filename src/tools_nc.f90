!----------------------------------------------------------------------
! Module: tools_nc
!> Purpose: NetCDF routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module tools_nc

use netcdf
use tools_display, only: msgerror
implicit none

! NetCDF type for floats
integer :: ncfloat !< NetCDF type for floats

private
public :: ncfloat
public :: ncfloat_init,ncerr

contains

!----------------------------------------------------------------------
! Subroutine: ncfloat_init
!> Purpose: initialize NetCDF type for floats
!----------------------------------------------------------------------
subroutine ncfloat_init

implicit none

! Check float size
if (kind(1.0)==4) then
   ncfloat = nf90_float
elseif (kind(1.0)==8) then
   ncfloat = nf90_double
else
   call msgerror('unknown real kind for NetCDF floats')
end if

end subroutine ncfloat_init

!----------------------------------------------------------------------
! Subroutine: ncerr
!> Purpose: handle NetCDF error
!----------------------------------------------------------------------
subroutine ncerr(subr,info)

implicit none

! Passed variables
character(len=*),intent(in) :: subr !< Calling subroutine
integer,intent(in) :: info          !< Info index

! Check status
if (info/=nf90_noerr) call msgerror('in '//trim(subr)//': '//trim(nf90_strerror(info)))

end subroutine ncerr

end module tools_nc

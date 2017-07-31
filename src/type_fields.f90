!----------------------------------------------------------------------
! Module: type_fields
!> Purpose: fields derived types
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module type_fields

implicit none

! Full field derived type
type fldtype
   real,allocatable :: val(:,:)  !< Global field data
   real,allocatable :: vala(:,:) !< Local field data, halo A
end type fldtype

! Reduced field derived type
type alphatype
   real,allocatable :: val(:)  !< Global subgrid variable data
   real,allocatable :: vala(:) !< Local subgrid variable data, halo A
   real,allocatable :: valb(:) !< Local subgrid variable data, halo B
   real,allocatable :: valc(:) !< Local subgrid variable data, halo C
end type alphatype

! Buffer derived type
type buftype
   real,allocatable :: val(:) !< Buffer data
end type buftype

private
public :: fldtype,alphatype,buftype

end module type_fields

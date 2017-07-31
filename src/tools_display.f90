!----------------------------------------------------------------------
! Module: tools_display
!> Purpose: display variables
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module tools_display

use type_mpl, only: mpl,mpl_abort
implicit none

! Display colors
character(len=1024) :: black     !< Black color code
character(len=1024) :: err       !< Error color code
character(len=1024) :: wng       !< Warning color code

! Progression display
integer :: ddis                  !< Progression display step

private
public :: black,err,wng,ddis
public :: msgerror,msgwarning,prog_init,prog_print

contains

!----------------------------------------------------------------------
! Subroutine: msgerror
!> Purpose: print error message and stop
!----------------------------------------------------------------------
subroutine msgerror(message)

implicit none

! Passed variables
character(len=*),intent(in) :: message !< Message

! Clean MPL abort
call mpl_abort(trim(err)//'!!! Error: '//trim(message)//trim(black))

end subroutine msgerror

!----------------------------------------------------------------------
! Subroutine: msgwarning
!> Purpose: print warning message
!----------------------------------------------------------------------
subroutine msgwarning(message)

implicit none

! Passed variables
character(len=*),intent(in) :: message !< Message

! Print warning message
write(mpl%unit,'(a)') trim(wng)//'!!! Warning: '//trim(message)//trim(black)

end subroutine msgwarning

!----------------------------------------------------------------------
! Subroutine: prog_init
!> Purpose: initialize progression display
!----------------------------------------------------------------------
subroutine prog_init(progint,done)

implicit none

! Passed variables
integer,intent(out) :: progint                    !< Progression integer
logical,dimension(:),intent(out),optional :: done !< Progression logical array

! Print message
write(mpl%unit,'(i3,a)',advance='no') 0,'%'

! Initialization
progint = ddis
if (present(done)) done = .false.

end subroutine prog_init

!----------------------------------------------------------------------
! Subroutine: prog_print
!> Purpose: print progression display
!----------------------------------------------------------------------
subroutine prog_print(progint,done)

implicit none

! Passed variables
integer,intent(inout) :: progint        !< Progression integer
logical,dimension(:),intent(in) :: done !< Progression logical array

! Local variables
real :: prog

! Print message
prog = 100.0*float(count(done))/float(size(done))
if (int(prog)>progint) then
   write(mpl%unit,'(i3,a)',advance='no') int(progint),'% '
   progint = progint+ddis
end if

end subroutine prog_print

end module tools_display

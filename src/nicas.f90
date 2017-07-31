!----------------------------------------------------------------------
! Program: nicas
!> Purpose: compute NICAS correlation model parameters
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
program nicas

use model_interface, only: model_coord
use module_mpi, only: compute_mpi
use module_namelist, only: model,new_param,new_mpi,check_adjoints,check_pos_def,check_mpi,check_dirac,namread,nproc,mpicom
use module_normalization, only: compute_normalization
use module_parameters, only: compute_parameters
use module_test, only: test_adjoints,test_pos_def,test_mpi,test_dirac
use tools_display, only: msgerror
use tools_nc, only: ncfloat_init
use tools_missing, only: isnotmsr,msr
use type_mpl, only: mpl,mpl_start,mpl_bcast,mpl_end
use type_sdata, only: sdatatype,sdata_read_param,sdata_read_local,sdata_read_mpi,sdata_write_param,sdata_write_mpi
use type_timer, only: timertype,timer_start,timer_intermediate,timer_end
implicit none

! Local variables
type(sdatatype) :: sdata
type(timertype) :: timer

!----------------------------------------------------------------------
! Initialize MPL
!----------------------------------------------------------------------

call mpl_start

!----------------------------------------------------------------------
! Header
!----------------------------------------------------------------------

write(mpl%unit,'(a)') '-------------------------------------------------------------------'
write(mpl%unit,'(a)') '--- You are running nicas -----------------------------------------'
write(mpl%unit,'(a)') '--- Author: Benjamin Menetrier ------------------------------------'
write(mpl%unit,'(a)') '--- Copyright © 2017 METEO-FRANCE------------------ ---------------'
write(mpl%unit,'(a)') '-------------------------------------------------------------------'
write(mpl%unit,'(a,i2,a,i2,a)') '--- Parallelization with ',mpl%nproc,' MPI tasks and ',mpl%nthread,' OpenMP threads'

!----------------------------------------------------------------------
! Timer
!----------------------------------------------------------------------

call timer_start(timer)

!----------------------------------------------------------------------
! Read and check namelist
!----------------------------------------------------------------------

write(mpl%unit,'(a)') '-------------------------------------------------------------------'
write(mpl%unit,'(a)') '--- Read and check namelist'

call namread

!----------------------------------------------------------------------
! Initialize random seeds and constants
!----------------------------------------------------------------------

write(mpl%unit,'(a)') '-------------------------------------------------------------------'
write(mpl%unit,'(a)') '--- Initialize constants'

call ncfloat_init

!----------------------------------------------------------------------
! Initialize OpenMP
!----------------------------------------------------------------------


!----------------------------------------------------------------------
! Initialize sampling
!----------------------------------------------------------------------

write(mpl%unit,'(a)') '-------------------------------------------------------------------'
write(mpl%unit,'(a,i5,a)') '--- Initialize sampling'

! Get coordinates
call model_coord(sdata)

! Initialize sdata parameters from namelist
sdata%nproc = max(nproc,1)
sdata%mpicom = mpicom

if (new_param) then
   !----------------------------------------------------------------------
   ! Compute NICAS parameters
   !----------------------------------------------------------------------

   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Compute NICAS parameters'

   call compute_parameters(sdata)

   !----------------------------------------------------------------------
   ! Compute NICAS normalization
   !----------------------------------------------------------------------

   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Compute NICAS normalization'

   call compute_normalization(sdata)

   !----------------------------------------------------------------------
   ! Write NICAS parameters
   !----------------------------------------------------------------------

   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Write NICAS parameters'

   call sdata_write_param(sdata)
else
   !----------------------------------------------------------------------
   ! Read NICAS parameters
   !----------------------------------------------------------------------

   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Read NICAS parameters'

   call sdata_read_param(sdata)
end if

if (new_mpi) then
   !----------------------------------------------------------------------
   ! Read NICAS local distribution
   !----------------------------------------------------------------------

   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Read NICAS local distribution'

   call sdata_read_local(sdata)

   !----------------------------------------------------------------------
   ! Compute NICAS MPI distribution
   !----------------------------------------------------------------------

   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Compute NICAS MPI distribution'

   sdata%mpicom = mpicom
   call compute_mpi(sdata)

   !----------------------------------------------------------------------
   ! Write NICAS MPI distribution
   !----------------------------------------------------------------------

   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Write NICAS MPI distribution'

   call sdata_write_mpi(sdata)
else
   !----------------------------------------------------------------------
   ! Read NICAS local distribution
   !----------------------------------------------------------------------

   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Read NICAS local distribution'

   call sdata_read_local(sdata)

   !----------------------------------------------------------------------
   ! Read NICAS MPI distribution
   !----------------------------------------------------------------------

   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Read NICAS MPI distribution'

   call sdata_read_mpi(sdata)
end if

if (check_adjoints) then
   !----------------------------------------------------------------------
   ! Test adjoints
   !----------------------------------------------------------------------

   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Test adjoints'

   call test_adjoints(sdata)
end if

if (check_pos_def) then
   !----------------------------------------------------------------------
   ! Test NICAS positive definiteness
   !----------------------------------------------------------------------

   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Test NICAS positive definiteness'

   call test_pos_def(sdata)
end if

if (check_mpi.and.(nproc>0)) then
   !----------------------------------------------------------------------
   ! Test single/multi-procs equivalence
   !----------------------------------------------------------------------

   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Test single/multi-procs equivalence'

  call test_mpi(sdata)
end if

if (check_dirac) then
   !----------------------------------------------------------------------
   ! Apply NICAS to diracs
   !----------------------------------------------------------------------

   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Apply NICAS to diracs'

   call test_dirac(sdata)
end if

!----------------------------------------------------------------------
! Execution stats
!----------------------------------------------------------------------

write(mpl%unit,'(a)') '-------------------------------------------------------------------'
write(mpl%unit,'(a)') '--- Execution stats'

call timer_end(timer)

write(mpl%unit,'(a)') '-------------------------------------------------------------------'

!----------------------------------------------------------------------
! Finalize MPL
!----------------------------------------------------------------------

call mpl_end

end program nicas

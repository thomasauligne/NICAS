!----------------------------------------------------------------------
! Module: module_parameters_interp.f90
!> Purpose: compute NICAS parameters (interpolation)
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module module_parameters_interp

use model_interface, only: model_read,model_write
use module_namelist, only: datadir,model,levs,mask_check,ntry,nrep,lsqrt,Lbh_file,Lbh,Lbv_file,Lbv,resol
use netcdf
use omp_lib
use tools_const, only: pi,req,deg2rad,rad2deg,sphere_dist,vector_product,vector_triple_product
use tools_display, only: msgerror,prog_init,prog_print
use tools_missing, only: msvali,msvalr,msi,msr,isnotmsr,isnotmsi
use tools_nc, only: ncfloat,ncerr
use type_ctree, only: ctreetype,create_ctree,find_nearest_neighbors,delete_ctree
use type_linop, only: linoptype,linop_alloc,linop_dealloc
use type_mpl, only: mpl,mpl_bcast,mpl_recv,mpl_send
use type_randgen, only: initialize_sampling,rand_integer
use type_sdata, only: sdatatype
implicit none

real,parameter :: S_inf = 1.0e-12 !< Minimum value for the interpolation coefficients

private
public :: compute_interp_h,compute_interp_v,compute_interp_s,esmf_grid

contains

!----------------------------------------------------------------------
! Subroutine: compute_interp_h
!> Purpose: compute basic horizontal interpolation
!----------------------------------------------------------------------
subroutine compute_interp_h(sdata)

implicit none

! Passed variables
type(sdatatype),intent(inout) :: sdata !< Sampling data

! Local variables
integer :: ncid_base,n_s_id,row_id,col_id,S_id
integer :: ic0,ic1,get_pid,i_s,ibnd,jc0,jc1,progint,il0i
integer :: iproc,i_s_s(mpl%nproc),i_s_e(mpl%nproc),n_s_loc(mpl%nproc),i_s_loc
integer,allocatable :: mask_ctree(:)
real :: dum(1)
real :: renorm(sdata%nc0)
real,allocatable :: x(:),y(:),z(:),v1(:),v2(:),va(:),vp(:),t(:)
logical,allocatable :: done(:),valid(:),validg(:,:),missing(:)
character(len=4) :: myprocchar
character(len=8) :: pidchar
character(len=1024) :: samplingname_base,interpname_base
character(len=1024) :: subr = 'compute_interp'
type(linoptype) :: hbase,htmp
type(ctreetype) :: ctree

! Get pid to build names (for the case where several nicas run simultaneously)
call fgetpid(get_pid)
write(pidchar,'(i8.8)') get_pid

! Build names
write(myprocchar,'(i4.4)') mpl%myproc
samplingname_base = 'sampling_data_base_'//pidchar//'_'//myprocchar//'.nc'
interpname_base = 'interp_data_base_'//pidchar//'_'//myprocchar//'.nc'

! Write grid data for ESMF_RegridWeightGen
call esmf_grid(sdata%nc1,sdata%lon(sdata%ic1_to_ic0(1:sdata%nc1)),sdata%lat(sdata%ic1_to_ic0(1:sdata%nc1)),''//trim(samplingname_base)//'')

! Compute interpolation weights with ESMF_RegridWeightGen
if (all(sdata%area<4.0*pi)) then
   ! LAM, lat/lon
   call execute_command_line('ESMF_RegridWeightGen --no_log -i -r -s '//trim(datadir)//'/'//trim(samplingname_base)//' --src_type ESMF -d '//trim(datadir)//'/grid_SCRIP.nc --dst_type SCRIP -m bilinear -w '//trim(datadir)//'/'//trim(interpname_base)//' > ESMF_log')
else
   ! Global model
   if (isnotmsi(sdata%nlon).and.isnotmsi(sdata%nlon)) then
      ! Lat/lon grid
      call execute_command_line('ESMF_RegridWeightGen --no_log -i -s '//trim(datadir)//'/'//trim(samplingname_base)//' --src_type ESMF -d '//trim(datadir)//'/grid_SCRIP.nc --dst_type SCRIP -m bilinear -w '//trim(datadir)//'/'//trim(interpname_base)//' > ESMF_log')
   else
      ! Unstructured grid
      call execute_command_line('ESMF_RegridWeightGen --no_log -i -s '//trim(datadir)//'/'//trim(samplingname_base)//' --src_type ESMF -d '//trim(datadir)//'/grid_SCRIP.nc --dst_type ESMF -m bilinear -w '//trim(datadir)//'/'//trim(interpname_base)//' > ESMF_log')
   end if
end if

! Read interpolation data size
call ncerr(subr,nf90_open(trim(datadir)//'/'//trim(interpname_base),nf90_nowrite,ncid_base))
call ncerr(subr,nf90_inq_dimid(ncid_base,'n_s',n_s_id))
call ncerr(subr,nf90_inquire_dimension(ncid_base,n_s_id,len=hbase%n_s))

! Allocation
call linop_alloc(hbase)


! Read interpolation data, close and delete files
call ncerr(subr,nf90_inq_varid(ncid_base,'row',row_id))
call ncerr(subr,nf90_inq_varid(ncid_base,'col',col_id))
call ncerr(subr,nf90_inq_varid(ncid_base,'S',S_id))
call ncerr(subr,nf90_get_var(ncid_base,row_id,hbase%row))
call ncerr(subr,nf90_get_var(ncid_base,col_id,hbase%col))
call ncerr(subr,nf90_get_var(ncid_base,S_id,hbase%S))
call ncerr(subr,nf90_close(ncid_base))
!call execute_command_line('rm -f '//trim(datadir)//'/'//trim(interpname_base)//' ESMF_log PET*.RegridWeightGen.Log')

! Conversion
hbase%row = sdata%ic0(hbase%row)

! Allocation
allocate(sdata%h(sdata%nl0i))

do il0i=1,sdata%nl0i
   ! Allocation and copy
   htmp%n_s = hbase%n_s
   call linop_alloc(htmp)
   htmp%row = hbase%row
   htmp%col = hbase%col
   htmp%S = hbase%S
   allocate(valid(htmp%n_s))

   ! Check interpolation coefficient
   valid = (htmp%S>S_inf)

   ! Check mask
   do i_s=1,htmp%n_s
      if (valid(i_s)) then
         ic0 = htmp%row(i_s)
         jc1 = htmp%col(i_s)
         jc0 = sdata%ic1_to_ic0(jc1)
         valid(i_s) = sdata%mask(ic0,il0i).and.sdata%mask(jc0,il0i)
      end if
   end do

   if (mask_check) then
      ! MPI splitting
      do iproc=1,mpl%nproc
         i_s_s(iproc) = (iproc-1)*(htmp%n_s/mpl%nproc+1)+1
         i_s_e(iproc) = min(iproc*(htmp%n_s/mpl%nproc+1),htmp%n_s)
         n_s_loc(iproc) = i_s_e(iproc)-i_s_s(iproc)+1
      end do

      ! Allocation
      allocate(done(n_s_loc(mpl%myproc)))

      ! Check that interpolations are not crossing mask boundaries
      write(mpl%unit,'(a10,a,i3,a)',advance='no') '','Level ',levs(il0i),': '
      call prog_init(progint,done)
      !$omp parallel do private(i_s_loc,i_s,x,y,z,v1,v2,va,vp,t,ic0,jc1,jc0)
      do i_s_loc=1,n_s_loc(mpl%myproc)
         ! Indices
         i_s = i_s_s(mpl%myproc)+i_s_loc-1

         if (valid(i_s)) then
            ! Allocation
            allocate(x(2))
            allocate(y(2))
            allocate(z(2))
            allocate(v1(3))
            allocate(v2(3))
            allocate(va(3))
            allocate(vp(3))
            allocate(t(4))

            ! Indices
            ic0 = htmp%row(i_s)
            jc1 = htmp%col(i_s)
            jc0 = sdata%ic1_to_ic0(jc1)

            ! Transform to cartesian coordinates
            call trans(2,sdata%lat((/ic0,jc0/)),sdata%lon((/ic0,jc0/)),x,y,z)

            ! Compute arc orthogonal vector
            v1 = (/x(1),y(1),z(1)/)
            v2 = (/x(2),y(2),z(2)/)
            call vector_product(v1,v2,va)

            ! Check if arc is crossing boundary arcs
            do ibnd=1,sdata%nbnd(il0i)
               call vector_product(va,sdata%vbnd(:,ibnd,il0i),vp)
               v1 = (/x(1),y(1),z(1)/)
               call vector_triple_product(v1,va,vp,t(1))
               v1 = (/x(2),y(2),z(2)/)
               call vector_triple_product(v1,va,vp,t(2))
               v1 = (/sdata%xbnd(1,ibnd,il0i),sdata%ybnd(1,ibnd,il0i),sdata%zbnd(1,ibnd,il0i)/)
               call vector_triple_product(v1,sdata%vbnd(:,ibnd,il0i),vp,t(3))
               v1 = (/sdata%xbnd(2,ibnd,il0i),sdata%ybnd(2,ibnd,il0i),sdata%zbnd(2,ibnd,il0i)/)
               call vector_triple_product(v1,sdata%vbnd(:,ibnd,il0i),vp,t(4))
               t(1) = -t(1)
               t(3) = -t(3)
               if (all(t>0).or.(all(t<0))) then
                  valid(i_s) = .false.
                  exit
               end if
            end do

            ! Memory release
            deallocate(x)
            deallocate(y)
            deallocate(z)
            deallocate(v1)
            deallocate(v2)
            deallocate(va)
            deallocate(vp)
            deallocate(t)
         end if

         ! Print progression
         done(i_s_loc) = .true.
         call prog_print(progint,done)
      end do
      !$omp end parallel do
      write(mpl%unit,'(a)') '100%'

      ! Communication
      if (mpl%main) then
         ! Allocation
         allocate(validg(maxval(n_s_loc),mpl%nproc))

         do iproc=1,mpl%nproc
            if (iproc==mpl%ioproc) then
               ! Copy data
               validg(1:n_s_loc(iproc),iproc) = valid(i_s_s(iproc):i_s_e(iproc))
            else
               ! Receive data on ioproc
               call mpl_recv(n_s_loc(iproc),validg(1:n_s_loc(iproc),iproc),iproc,mpl%tag)
            end if
         end do

         ! Format data
         do iproc=1,mpl%nproc
            valid(i_s_s(iproc):i_s_e(iproc)) = validg(1:n_s_loc(iproc),iproc)
         end do

         ! Release memory
         deallocate(validg)
      else
         ! Send data to ioproc
         call mpl_send(n_s_loc(mpl%myproc),valid(i_s_s(mpl%myproc):i_s_e(mpl%myproc)),mpl%ioproc,mpl%tag)
      end if
      mpl%tag = mpl%tag+1

      ! Broadcast
      call mpl_bcast(valid,mpl%ioproc)

      ! Release memory
      deallocate(done)
   else
      write(mpl%unit,'(a10,a,i3)') '','Level ',levs(il0i)
   end if

   ! Renormalization
   renorm = 0.0
   do i_s=1,htmp%n_s
      if (valid(i_s)) renorm(htmp%row(i_s)) = renorm(htmp%row(i_s))+htmp%S(i_s)
   end do

   ! Initialize object
   sdata%h(il0i)%prefix = 'h'
   sdata%h(il0i)%n_src = sdata%nc1
   sdata%h(il0i)%n_dst = sdata%nc0
   sdata%h(il0i)%n_s = count(valid)
   call linop_alloc(sdata%h(il0i))
   sdata%h(il0i)%n_s = 0
   do i_s=1,htmp%n_s
      if (valid(i_s)) then
         sdata%h(il0i)%n_s = sdata%h(il0i)%n_s+1
         sdata%h(il0i)%row(sdata%h(il0i)%n_s) = htmp%row(i_s)
         sdata%h(il0i)%col(sdata%h(il0i)%n_s) = htmp%col(i_s)
         sdata%h(il0i)%S(sdata%h(il0i)%n_s) = htmp%S(i_s)/renorm(htmp%row(i_s))
      end if
   end do

   ! Release memory
   call linop_dealloc(htmp)
   deallocate(valid)

   ! Allocation
   allocate(missing(sdata%nc0))

   ! Count points that are not interpolated
   missing = .false.
   do ic0=1,sdata%nc0
      if (sdata%mask(ic0,il0i)) missing(ic0) = .true.
   end do
   do i_s=1,sdata%h(il0i)%n_s
      missing(sdata%h(il0i)%row(i_s)) = .false.
   end do
   if (count(missing)>0) then
      ! Allocate temporary arrays
      htmp%n_s = sdata%h(il0i)%n_s
      call linop_alloc(htmp)

      ! Fill temporary arrays
      htmp%row = sdata%h(il0i)%row
      htmp%col = sdata%h(il0i)%col
      htmp%S = sdata%h(il0i)%S

      ! Reallocate interpolation
      call linop_dealloc(sdata%h(il0i))
      sdata%h(il0i)%n_s = sdata%h(il0i)%n_s+count(missing)
      call linop_alloc(sdata%h(il0i))

      ! Fill permanent arrays
      sdata%h(il0i)%row(1:htmp%n_s) = htmp%row
      sdata%h(il0i)%col(1:htmp%n_s) = htmp%col
      sdata%h(il0i)%S(1:htmp%n_s) = htmp%S

      ! Compute cover tree
      allocate(mask_ctree(sdata%nc1))
      mask_ctree = 0
      do ic1=1,sdata%nc1
         if (sdata%mask(sdata%ic1_to_ic0(ic1),il0i)) mask_ctree(ic1) = 1
      end do
      ctree = create_ctree(sdata%nc1,dble(sdata%lon(sdata%ic1_to_ic0)),dble(sdata%lat(sdata%ic1_to_ic0)),mask_ctree)
      deallocate(mask_ctree)

      ! Compute nearest neighbors
      do ic0=1,sdata%nc0
         if (missing(ic0)) then
            htmp%n_s = htmp%n_s+1
            sdata%h(il0i)%row(htmp%n_s) = ic0
            call find_nearest_neighbors(ctree,dble(sdata%lon(ic0)),dble(sdata%lat(ic0)),1,sdata%h(il0i)%col(htmp%n_s:htmp%n_s),dum)
            sdata%h(il0i)%S(htmp%n_s) = 1.0
         end if
      end do

      ! Release memory
      call linop_dealloc(htmp)
      call delete_ctree(ctree)
   end if

   ! Release memory
   deallocate(missing)
end do

! Release memory
call linop_dealloc(hbase)

end subroutine compute_interp_h

!----------------------------------------------------------------------
! Subroutine: compute_interp_v
!> Purpose: compute vertical interpolation
!----------------------------------------------------------------------
subroutine compute_interp_v(sdata)

implicit none

! Passed variables
type(sdatatype),intent(inout) :: sdata !< Sampling data

! Local variables
integer :: il0,jl0,il1,il0inf,il0sup,il0i,i_s,ic0,ic1
logical,allocatable :: valid(:)
type(linoptype) :: vbase,vtmp

! Linear interpolation
vbase%n_s = sdata%nl1
il0inf = 1
do il0=1,sdata%nl0
   if (sdata%llev(il0)) then
      il0sup = il0
      do jl0=il0inf+1,il0sup-1
         vbase%n_s = vbase%n_s+2
      end do
      il0inf = il0
   end if
end do
call linop_alloc(vbase)
do il1=1,sdata%nl1
   il0 = sdata%il1_to_il0(il1)
   vbase%row(il1) = il0
   vbase%col(il1) = il0
   vbase%S(il1) = 1.0
end do
vbase%n_s = sdata%nl1
il0inf = 1
do il0=1,sdata%nl0
   if (sdata%llev(il0)) then
      il0sup = il0
      do jl0=il0inf+1,il0sup-1
         vbase%n_s = vbase%n_s+1
         vbase%row(vbase%n_s) = jl0
         vbase%col(vbase%n_s) = il0inf
         vbase%S(vbase%n_s) = abs(sdata%vunit(il0sup)-sdata%vunit(jl0))/abs(sdata%vunit(il0sup)-sdata%vunit(il0inf))

         vbase%n_s = vbase%n_s+1
         vbase%row(vbase%n_s) = jl0
         vbase%col(vbase%n_s) = il0sup
         vbase%S(vbase%n_s) = abs(sdata%vunit(jl0)-sdata%vunit(il0inf))/abs(sdata%vunit(il0sup)-sdata%vunit(il0inf))
      end do
      il0inf = il0
   end if
end do

! Allocation
allocate(sdata%v(sdata%nl0i))
allocate(valid(vbase%n_s))

do il0i=1,sdata%nl0i
   ! Initialize vertical interpolation
   sdata%v(il0i)%prefix = 'v'
   sdata%v(il0i)%n_src = sdata%nl1
   sdata%v(il0i)%n_dst = sdata%nl0

   if (sdata%nl0i==1) then
      ! Copy basic vertical interpolation
      sdata%v(il0i)%n_s = vbase%n_s
      call linop_alloc(sdata%v(il0i))
      sdata%v(il0i)%row = vbase%row
      sdata%v(il0i)%col = vbase%col
      sdata%v(il0i)%S = vbase%S
   else
      ! Check valid operations
      do i_s=1,vbase%n_s
         valid(i_s) = (vbase%row(i_s)<=il0i).and.(vbase%col(i_s)<=il0i)
      end do

      ! Copy valid operations
      vtmp%n_s = count(valid)
      call linop_alloc(vtmp)
      vtmp%n_s = 0
      do i_s=1,vbase%n_s
         if (valid(i_s)) then
            vtmp%n_s = vtmp%n_s+1
            vtmp%row(vtmp%n_s) = vbase%row(i_s)
            vtmp%col(vtmp%n_s) = vbase%col(i_s)
            vtmp%S(vtmp%n_s) = vbase%S(i_s)
         end if
      end do

      ! Add missing levels to the interpolation
      sdata%v(il0i)%n_s = vtmp%n_s+(il0i-maxval(vtmp%col))
      call linop_alloc(sdata%v(il0i))
      sdata%v(il0i)%row(1:vtmp%n_s) = vtmp%row
      sdata%v(il0i)%col(1:vtmp%n_s) = vtmp%col
      sdata%v(il0i)%S(1:vtmp%n_s) = vtmp%S
      do jl0=maxval(vtmp%col)+1,il0i
         vtmp%n_s = vtmp%n_s+1
         sdata%v(il0i)%row(vtmp%n_s) = jl0
         sdata%v(il0i)%col(vtmp%n_s) = jl0
         sdata%v(il0i)%S(vtmp%n_s) = 1.0
      end do

      ! Release memory
      call linop_dealloc(vtmp)
   end if

   ! Conversion
   sdata%v(il0i)%col = sdata%il0_to_il1(sdata%v(il0i)%col)
end do

! Release memory
call linop_dealloc(vbase)
deallocate(valid)

! Find the bottom for each point of S1
allocate(sdata%vbot(sdata%nc1))
!$omp parallel do private(ic1,ic0,il0)
do ic1=1,sdata%nc1
   ic0 = sdata%ic1_to_ic0(ic1)
   il0 = 1
   do while (sdata%mask(ic0,il0).and.(il0<sdata%nl0))
      il0 = il0+1
   end do
   sdata%vbot(ic1) = min(il0,sdata%nl0i)
end do
!$omp end parallel do

end subroutine compute_interp_v

!----------------------------------------------------------------------
! Subroutine: compute_interp_s
!> Purpose: compute horizontal subsampling interpolation
!----------------------------------------------------------------------
subroutine compute_interp_s(sdata)

implicit none

! Passed variables
type(sdatatype),intent(inout) :: sdata !< Sampling data

! Local variables
integer :: ncid,n_s_id,row_id,col_id,S_id
integer :: ic1,ic2,il1,get_pid,i_s,ibnd,il0,ic0,jc2,jc1,jc0,progint
integer :: iproc,i_s_s(mpl%nproc),i_s_e(mpl%nproc),n_s_loc(mpl%nproc),i_s_loc
integer,allocatable :: mask_ctree(:)
real :: dum(1)
real :: renorm(sdata%nc1)
real,allocatable :: x(:),y(:),z(:),v1(:),v2(:),va(:),vp(:),t(:)
logical,allocatable :: done(:),valid(:),validg(:,:),missing(:)
character(len=3) :: ilschar
character(len=4) :: myprocchar
character(len=8) :: pidchar
character(len=1024) :: samplingname_base,samplingname,interpname
character(len=1024) :: subr = 'compute_interp'
type(linoptype) :: stmp
type(ctreetype) :: ctree

! Allocation
allocate(sdata%s(sdata%nl1))

! Get pid to build names (for the case where several nicas run simultaneously)
call fgetpid(get_pid)
write(pidchar,'(i8.8)') get_pid

! Build names
write(myprocchar,'(i4.4)') mpl%myproc
samplingname_base = 'sampling_data_base_'//pidchar//'_'//myprocchar//'.nc'

do il1=1,sdata%nl1
   ! Build names
   write(ilschar,'(i3.3)') il1
   samplingname = 'sampling_data_'//ilschar//'_'//pidchar//'_'//myprocchar//'.nc'
   interpname = 'interp_data_'//ilschar//'_'//pidchar//'_'//myprocchar//'.nc'

   ! Write grid data for ESMF_RegridWeightGen
   call esmf_grid(sdata%nc2(il1),sdata%lon(sdata%ic2il1_to_ic0(1:sdata%nc2(il1),il1)),sdata%lat(sdata%ic2il1_to_ic0(1:sdata%nc2(il1),il1)),''//trim(samplingname)//'')

   ! Compute interpolation weights with ESMF_RegridWeightGen
   if (all(sdata%area<4.0*pi)) then
      ! LAM
      call execute_command_line('ESMF_RegridWeightGen --no_log -i -r -s '//trim(datadir)//'/'//trim(samplingname)//' --src_type ESMF -d '//trim(datadir)//'/'//trim(samplingname_base)//' --dst_type ESMF -m bilinear -w '//trim(datadir)//'/'//trim(interpname)//' > ESMF_log')
   else
      ! Global model
      call execute_command_line('ESMF_RegridWeightGen --no_log -s '//trim(datadir)//'/'//trim(samplingname)//' --src_type ESMF -d '//trim(datadir)//'/'//trim(samplingname_base)//' --dst_type ESMF -m bilinear -w '//trim(datadir)//'/'//trim(interpname)//' > ESMF_log')
   end if

   ! Read interpolation data size
   call ncerr(subr,nf90_open(trim(datadir)//'/'//trim(interpname),nf90_nowrite,ncid))
   call ncerr(subr,nf90_inq_dimid(ncid,'n_s',n_s_id))
   call ncerr(subr,nf90_inquire_dimension(ncid,n_s_id,len=stmp%n_s))

   ! Allocation
   call linop_alloc(stmp)
   allocate(valid(stmp%n_s))

   ! Read interpolation data, close and delete files
   call ncerr(subr,nf90_inq_varid(ncid,'row',row_id))
   call ncerr(subr,nf90_inq_varid(ncid,'col',col_id))
   call ncerr(subr,nf90_inq_varid(ncid,'S',S_id))
   call ncerr(subr,nf90_get_var(ncid,row_id,stmp%row))
   call ncerr(subr,nf90_get_var(ncid,col_id,stmp%col))
   call ncerr(subr,nf90_get_var(ncid,S_id,stmp%S))
   call ncerr(subr,nf90_close(ncid))
!   call execute_command_line('rm -f '//trim(datadir)//'/'//trim(samplingname)//' '//trim(datadir)//'/'//trim(interpname)//' ESMF_log PET*.RegridWeightGen.Log')

   ! Check interpolation coefficient
   valid = .not.(stmp%S<S_inf)

   if (mask_check) then
      ! MPI splitting
      do iproc=1,mpl%nproc
         i_s_s(iproc) = (iproc-1)*(stmp%n_s/mpl%nproc+1)+1
         i_s_e(iproc) = min(iproc*(stmp%n_s/mpl%nproc+1),stmp%n_s)
         n_s_loc(iproc) = i_s_e(iproc)-i_s_s(iproc)+1
      end do

      ! Allocation
      allocate(done(n_s_loc(mpl%myproc)))

      ! Check that interpolations are not crossing mask boundaries
      write(mpl%unit,'(a10,a,i3,a)',advance='no') '','Level ',levs(sdata%il1_to_il0(il1)),': '
      call prog_init(progint,done)
      !$omp parallel do private(i_s_loc,i_s,x,y,z,v1,v2,va,vp,t,ic1,ic0,jc2,jc1,jc0,il0)
      do i_s_loc=1,n_s_loc(mpl%myproc)
         ! Indices
         i_s = i_s_s(mpl%myproc)+i_s_loc-1

         if (valid(i_s)) then
            ! Allocation
            allocate(x(2))
            allocate(y(2))
            allocate(z(2))
            allocate(v1(3))
            allocate(v2(3))
            allocate(va(3))
            allocate(vp(3))
            allocate(t(4))

            ! Indices
            ic1 = stmp%row(i_s)
            ic0 = sdata%ic1_to_ic0(ic1)
            jc2 = stmp%col(i_s)
            jc1 = sdata%ic2il1_to_ic1(jc2,il1)
            jc0 = sdata%ic1_to_ic0(jc1)
            il0 = sdata%il1_to_il0(il1)

            ! Transform to cartesian coordinates
            call trans(2,sdata%lat((/ic0,jc0/)),sdata%lon((/ic0,jc0/)),x,y,z)

            ! Compute arc orthogonal vector
            v1 = (/x(1),y(1),z(1)/)
            v2 = (/x(2),y(2),z(2)/)
            call vector_product(v1,v2,va)

            ! Check if arc is crossing boundary arcs
            do ibnd=1,sdata%nbnd(il0)
               call vector_product(va,sdata%vbnd(:,ibnd,il0),vp)
               v1 = (/x(1),y(1),z(1)/)
               call vector_triple_product(v1,va,vp,t(1))
               v1 = (/x(2),y(2),z(2)/)
               call vector_triple_product(v1,va,vp,t(2))
               v1 = (/sdata%xbnd(1,ibnd,il0),sdata%ybnd(1,ibnd,il0),sdata%zbnd(1,ibnd,il0)/)
               call vector_triple_product(v1,sdata%vbnd(:,ibnd,il0),vp,t(3))
               v1 = (/sdata%xbnd(2,ibnd,il0),sdata%ybnd(2,ibnd,il0),sdata%zbnd(2,ibnd,il0)/)
               call vector_triple_product(v1,sdata%vbnd(:,ibnd,il0),vp,t(4))
               t(1) = -t(1)
               t(3) = -t(3)
               if (all(t>0).or.(all(t<0))) then
                  valid(i_s) = .false.
                  exit
               end if
            end do

            ! Memory release
            deallocate(x)
            deallocate(y)
            deallocate(z)
            deallocate(v1)
            deallocate(v2)
            deallocate(va)
            deallocate(vp)
            deallocate(t)
         end if

         ! Print progression
         done(i_s_loc) = .true.
         call prog_print(progint,done)
      end do
      !$omp end parallel do
      write(mpl%unit,'(a)') '100%'

      ! Communication
      if (mpl%main) then
         ! Allocation
         allocate(validg(maxval(n_s_loc),mpl%nproc))

         do iproc=1,mpl%nproc
            if (iproc==mpl%ioproc) then
               ! Copy data
               validg(1:n_s_loc(iproc),iproc) = valid(i_s_s(iproc):i_s_e(iproc))
            else
               ! Receive data on ioproc
               call mpl_recv(n_s_loc(iproc),validg(1:n_s_loc(iproc),iproc),iproc,mpl%tag)
            end if
         end do

         ! Format data
         do iproc=1,mpl%nproc
            valid(i_s_s(iproc):i_s_e(iproc)) = validg(1:n_s_loc(iproc),iproc)
         end do

         ! Release memory
         deallocate(validg)
      else
         ! Send data to ioproc
         call mpl_send(n_s_loc(mpl%myproc),valid(i_s_s(mpl%myproc):i_s_e(mpl%myproc)),mpl%ioproc,mpl%tag)
      end if
      mpl%tag = mpl%tag+1

      ! Broadcast
      call mpl_bcast(valid,mpl%ioproc)

      ! Release memory
      deallocate(done)
   else
      write(mpl%unit,'(a10,a,i3)') '','Level ',levs(sdata%il1_to_il0(il1))
   end if

   ! Renormalization
   renorm = 0.0
   do i_s=1,stmp%n_s
      if (valid(i_s)) renorm(stmp%row(i_s)) = renorm(stmp%row(i_s))+stmp%S(i_s)
   end do

   ! Initialize object
   sdata%s(il1)%prefix = 's'
   sdata%s(il1)%n_src = sdata%nc2(il1)
   sdata%s(il1)%n_dst = sdata%nc1
   sdata%s(il1)%n_s = count(valid)
   call linop_alloc(sdata%s(il1))
   sdata%s(il1)%n_s = 0
   do i_s=1,stmp%n_s
      if (valid(i_s)) then
         sdata%s(il1)%n_s = sdata%s(il1)%n_s+1
         sdata%s(il1)%row(sdata%s(il1)%n_s) = stmp%row(i_s)
         sdata%s(il1)%col(sdata%s(il1)%n_s) = stmp%col(i_s)
         sdata%s(il1)%S(sdata%s(il1)%n_s) = stmp%S(i_s)/renorm(stmp%row(i_s))
      end if
   end do

   ! Release memory
   call linop_dealloc(stmp)
   deallocate(valid)

   ! Allocation
   allocate(missing(sdata%nc1))

   ! Count points that are not interpolated
   missing = .false.
   do ic1=1,sdata%nc1
      if (sdata%mask(sdata%ic1_to_ic0(ic1),sdata%il1_to_il0(il1))) missing(ic1) = .true.
   end do
   do i_s=1,sdata%s(il1)%n_s
      missing(sdata%s(il1)%row(i_s)) = .false.
   end do
   if (count(missing)>0) then
      ! Allocate temporary arrays
      stmp%n_s = sdata%s(il1)%n_s
      call linop_alloc(stmp)

      ! Fill temporary arrays
      stmp%row = sdata%s(il1)%row
      stmp%col = sdata%s(il1)%col
      stmp%S = sdata%s(il1)%S

      ! Reallocate permanent arrays
      call linop_dealloc(sdata%s(il1))
      sdata%s(il1)%n_s = sdata%s(il1)%n_s+count(missing)
      call linop_alloc(sdata%s(il1))

      ! Fill permanent arrays
      sdata%s(il1)%row(1:stmp%n_s) = stmp%row
      sdata%s(il1)%col(1:stmp%n_s) = stmp%col
      sdata%s(il1)%S(1:stmp%n_s) = stmp%S

      ! Compute cover tree
      allocate(mask_ctree(sdata%nc2(il1)))
      mask_ctree = 0
      do ic2=1,sdata%nc2(il1)
         if (sdata%mask(sdata%ic1_to_ic0(sdata%ic2il1_to_ic1(ic2,il1)),sdata%il1_to_il0(il1))) mask_ctree(ic2) = 1
      end do
      ctree = create_ctree(sdata%nc2(il1),dble(sdata%lon(sdata%ic2il1_to_ic0(1:sdata%nc2(il1),il1))),dble(sdata%lat(sdata%ic2il1_to_ic0(1:sdata%nc2(il1),il1))),mask_ctree)
      deallocate(mask_ctree)

      ! Compute nearest neighbors
      do ic1=1,sdata%nc1
         if (missing(ic1)) then
            stmp%n_s = stmp%n_s+1
            sdata%s(il1)%row(stmp%n_s) = ic1
            call find_nearest_neighbors(ctree,dble(sdata%lon(sdata%ic1_to_ic0(ic1))),dble(sdata%lat(sdata%ic1_to_ic0(ic1))),1,sdata%s(il1)%col(stmp%n_s:stmp%n_s),dum)
            sdata%s(il1)%S(stmp%n_s) = 1.0
         end if
      end do

      ! Release memory
      call linop_dealloc(stmp)
      call delete_ctree(ctree)
   end if

   ! Release memory
   deallocate(missing)
end do

!call execute_command_line('rm -f '//trim(datadir)//'/'//trim(samplingname_base)//' ESMF_log PET*.RegridWeightGen.Log')

end subroutine compute_interp_s

!----------------------------------------------------------------------
! Subroutine: esmf_grid
!> Purpose: write ESMF grid in NetCDF format
!----------------------------------------------------------------------
subroutine esmf_grid(n,lon,lat,gridname)

implicit none

! Passed variables
integer,intent(in) :: n                 !< Sampling data
real,intent(in) :: lon(n)               !< Points longitudes
real,intent(in) :: lat(n)               !< Points latitudes
character(len=*),intent(in) :: gridname !< Grid file name

! Local variables
integer :: list(6*(n-2)),lptr(6*(n-2)),lend(n),lnew,near(n),next(n),info,nt,ltri(6,2*(n-2)),it
integer :: elementConn(3,2*(n-2)),numElementConn(2*(n-2)),elementMask(2*(n-2))
integer :: ncid,nodeCount_id,elementCount_id,maxNodePElement_id,coordDim_id,nodeCoords_id,elementConn_id,numElementConn_id,centerCoords_id,elementArea_id,elementMask_id
real :: x(n),y(n),z(n),dist(n),xc,yc,zc,normc,latc,lonc,rad,v1(3),v2(3),v3(3),areas
real :: nodeCoords(2,n),centerCoords(2,2*(n-2)),elementArea(2*(n-2))
character(len=1024) :: subr = 'esmf_grid'

! Transform to cartesian coordinates
call trans(n,lat,lon,x,y,z)

! Create mesh
list = 0
call trmesh(n,x,y,z,list,lptr,lend,lnew,near,next,dist,info)

! Create triangles list
call trlist(n,list,lptr,lend,6,nt,ltri,info)

! Fill ESMF mesh data
nodeCoords(1,:) = lon
nodeCoords(2,:) = lat
elementConn(:,1:nt) = ltri(1:3,1:nt)
numElementConn = 3
do it=1,nt
   xc = sum(x(ltri(1:3,it)))
   yc = sum(y(ltri(1:3,it)))
   zc = sum(z(ltri(1:3,it)))
   normc = sqrt(xc**2+yc**2+zc**2)
   xc = xc/normc
   yc = yc/normc
   zc = zc/normc
   call scoord(xc,yc,zc,latc,lonc,rad)
   centerCoords(1,it) = lonc
   centerCoords(2,it) = latc
   v1 = (/x(ltri(1,it)),y(ltri(1,it)),z(ltri(1,it))/)
   v2 = (/x(ltri(2,it)),y(ltri(2,it)),z(ltri(2,it))/)
   v3 = (/x(ltri(3,it)),y(ltri(3,it)),z(ltri(3,it))/)
   elementArea(it) = areas(v1,v2,v3)
end do
elementMask = 1

! Write ESMF data file
call ncerr(subr,nf90_create(trim(datadir)//'/'//trim(gridname),or(nf90_clobber,nf90_64bit_offset),ncid))
call ncerr(subr,nf90_def_dim(ncid,'nodeCount',n,nodeCount_id))
call ncerr(subr,nf90_def_dim(ncid,'elementCount',nt,elementCount_id))
call ncerr(subr,nf90_def_dim(ncid,'maxNodePElement',3,maxNodePElement_id))
call ncerr(subr,nf90_def_dim(ncid,'coordDim',2,coordDim_id))
call ncerr(subr,nf90_def_var(ncid,'nodeCoords',nf90_double,(/coordDim_id,nodeCount_id/),nodeCoords_id))
call ncerr(subr,nf90_put_att(ncid,nodeCoords_id,'units','radians'))
call ncerr(subr,nf90_def_var(ncid,'elementConn',nf90_int,(/maxNodePElement_id,elementCount_id/),elementConn_id))
call ncerr(subr,nf90_put_att(ncid,elementConn_id,'long_name','Node Indices that define the element / connectivity'))
call ncerr(subr,nf90_put_att(ncid,elementConn_id,'_FillValue',-1))
call ncerr(subr,nf90_def_var(ncid,'numElementConn',nf90_byte,(/elementCount_id/),numElementConn_id))
call ncerr(subr,nf90_put_att(ncid,numElementConn_id,'long_name','Number of nodes per element'))
call ncerr(subr,nf90_def_var(ncid,'centerCoords',nf90_double,(/coordDim_id,elementCount_id/),centerCoords_id))
call ncerr(subr,nf90_put_att(ncid,centerCoords_id,'units','radians'))
call ncerr(subr,nf90_def_var(ncid,'elementArea',nf90_double,(/elementCount_id/),elementArea_id))
call ncerr(subr,nf90_put_att(ncid,elementArea_id,'units','radians^2'))
call ncerr(subr,nf90_put_att(ncid,elementArea_id,'long_name','area weights'))
call ncerr(subr,nf90_def_var(ncid,'elementMask',nf90_int,(/elementCount_id/),elementMask_id))
call ncerr(subr,nf90_put_att(ncid,elementMask_id,'_FillValue',-9999))
call ncerr(subr,nf90_put_att(ncid,elementArea_id,'units','radians^2'))
call ncerr(subr,nf90_put_att(ncid,nf90_global,'gridType','unstructured'))
call ncerr(subr,nf90_put_att(ncid,nf90_global,'version','0.9'))
call ncerr(subr,nf90_enddef(ncid))
call ncerr(subr,nf90_put_var(ncid,nodeCoords_id,nodeCoords))
call ncerr(subr,nf90_put_var(ncid,elementConn_id,elementConn(:,1:nt)))
call ncerr(subr,nf90_put_var(ncid,numElementConn_id,numElementConn(1:nt)))
call ncerr(subr,nf90_put_var(ncid,centerCoords_id,centerCoords(:,1:nt)))
call ncerr(subr,nf90_put_var(ncid,elementArea_id,elementArea(1:nt)))
call ncerr(subr,nf90_put_var(ncid,elementMask_id,elementMask(1:nt)))
call ncerr(subr,nf90_close(ncid))

end subroutine esmf_grid

end module module_parameters_interp

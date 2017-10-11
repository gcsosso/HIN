program hin_dynamics

use, intrinsic :: iso_c_binding, only: C_NULL_CHAR, C_PTR, c_f_pointer
use xtc
use MOD_read_input 

implicit none

integer, parameter :: cart=3
integer :: i, j, t
integer :: eflag, fframe, stride, lframe, ns, it0, nat, resnum, idx, nsl
integer :: ntel, t0, isamp, nsamp, nframes, tmax, tt0, t0max, n_zmesh
integer :: NATOMS, STEP, STAT, counter, inif, delt, tl_ini, n_linear
integer, allocatable :: n_ws(:), list_ws(:,:)
integer, allocatable :: ntime(:), time0(:), normz(:)
real :: zmin, zmax, dz, dt, kvec, ttau, dummyp, kt(cart-1), dtime, box_trans(cart,cart), ktd(cart)
real :: prec, time, linear, tmin, ttau_min, ttau_min_TIME, zstep
real, allocatable :: pos(:,:), pos0(:,:), r2t(:,:), r2t_sisf(:,:), zmesh(:), sumx(:), sumy(:), sumxsq(:), sumxy(:)
real, allocatable :: r2t_TMP(:), r2t_sisf_TMP(:), bound(:,:), chi4(:,:), r2t_sisf_chi4(:,:), dc(:), relt(:)
real, parameter :: dcoeff=(1.0/(4.0*100.0))
character*3 :: hw_ex, switch_diff, switch_sisf, switch_chi4, td
character*5 :: resname
character*3, allocatable :: ws(:), sym(:)
character*100 :: sfile, tfile
type(C_PTR) :: xd_c
type(xdrfile), pointer :: xd

call read_input(eflag,sfile,tfile,fframe,stride,lframe,hw_ex,ns,ws,n_ws,zmin,zmax,dz, &
                it0,dt,kvec,ttau,switch_diff,switch_sisf,switch_chi4,linear,nsl,td)

! At the moment, just one species is allowed...
if (ns.gt.1) then
   write(*,*) "At the moment, just one species is allowed..."
   stop
endif   

! Get the components of the 2D q-vector (isotropic case)...
kt(:)=kvec/sqrt(2.0)

! If 3D (also isotropic case...)
ktd(:)=kvec/sqrt(3.0)

call read_gro(sfile,nat,sym,list_ws,n_ws,hw_ex,ns,resnum,resname,idx,dummyp,ws)

! Time & Z- stuff...

! Variables...
nframes=lframe-fframe
inif=fframe
nsamp=stride
dtime=dt*nsamp
tmax=int(real(nframes/2.0d0)/real(nsamp))
t0max=int(real(tmax)/real(it0))
n_zmesh=nsl !int((zmax-zmin)/dz)
zstep=(zmax-zmin)/real(nsl)

! Allocations...
allocate(r2t(tmax,n_zmesh),r2t_sisf(tmax,n_zmesh),chi4(tmax,n_zmesh),r2t_sisf_chi4(tmax,n_zmesh))
allocate(r2t_TMP(n_zmesh),r2t_sisf_TMP(n_zmesh))
allocate(ntime(tmax-1),time0(t0max))
allocate(zmesh(n_zmesh),bound(n_zmesh,2),normz(n_zmesh),dc(n_zmesh),relt(n_zmesh))
allocate(sumx(n_zmesh),sumy(n_zmesh),sumxsq(n_zmesh),sumxy(n_zmesh))

! Z-boundaries...
do i=1,n_zmesh
   j=i-1
   bound(i,1)=(real(j)*zstep)+zmin-(dz/2.0)
   bound(i,2)=bound(i,1)+dz
enddo

! Initializations...
ntel=0
t0=0
isamp=0
ntime(:)=0
r2t(:,:)=0.0
r2t_sisf(:,:)=0.0
chi4(:,:)=0.0
r2t_sisf_chi4(:,:)=0.0
sumx(:)=0.0
sumy(:)=0.0
sumxsq(:)=0.0
sumxy(:)=0.0

call read_first_xtc(tfile,STAT,NATOMS,nat,xd_c,xd,STEP,time,box_trans,pos,prec,cart,t0max,pos0)

! Read the whole thing
counter=0

do while ( STAT==0 )
    if (counter.gt.(inif-1)) then ! ignore some initial frames
       write(*,'(a,E12.6,a,i0,a,i0)') " Time (ps): ", time, "  Step: ", STEP, " Frame: ", counter
       if (mod(isamp,nsamp).eq.0) then
          ntel=ntel+1
          if (mod(ntel-1,it0).eq.0) then
             t0=t0+1
             if (t0.le.t0max) then
                tt0=t0
                time0(tt0)=ntel
                do i=1,nat
                   if (sym(i).eq.ws(1)) pos0(:,i+nat*(tt0-1))=pos(:,i)
                enddo
             endif
           endif
           do t=1,min(t0,t0max)
              delt=ntel-time0(t)+1
              if (delt.lt.tmax)then
                 ntime(delt)=ntime(delt)+1
                 ! Initialize stuff before the loop over the atoms
                 normz(:)=0
                 r2t_TMP(:)=0.0
                 r2t_sisf_TMP(:)=0.0
                 do i=1,nat
                    if (sym(i).eq.ws(1)) then
                       ! 2D
                       if (td.ne.'yes') then
                          do j=1,n_zmesh
                             if (pos(cart,i).ge.bound(j,1).and.pos(cart,i).lt.bound(j,2).and. & 
                                 pos0(cart,i+nat*(t-1)).ge.bound(j,1).and.pos0(cart,i+nat*(t-1)).lt.bound(j,2)) then
                                 ! Save stuff into tmp arrays running over the z-mesh...
                                 ! MSD
                                 r2t_TMP(j)=r2t_TMP(j)+(pos(1,i)-pos0(1,i+nat*(t-1)))*(pos(1,i)-pos0(1,i+nat*(t-1))) ! D_xx only ! 
                                 ! Remember: to get D, you should divide by 2 instead of 6 (! this is basically 1D, so /2*dim with dim=1) 
                                 ! SISF
                                 r2t_sisf_TMP(j)=r2t_sisf_TMP(j)+cos((pos(1,i)- &
                                                 pos0(1,i+nat*(t-1)))*kt(1)+(pos(2,i)-pos0(2,i+nat*(t-1)))*kt(2)) ! 2D !
                                 ! Norm: for every time window we have in principle
                                 ! a different number of atoms in the j-th z slice...
                                 normz(j)=normz(j)+1
                             endif
                          enddo
                       ! end 2D
                       else 
                       ! 3D
                          if (n_zmesh.ne.1) then
                             write(*,*) "You aske for a 3D calculation, but you're still doing slices!!"
                             stop
                          endif
                          do j=1,n_zmesh
                             ! 3D MSD
                             r2t_TMP(j)=r2t_TMP(j)+ dot_product( (pos(:,i)-pos0(:,i+nat*(t-1))) , (pos(:,i)-pos0(:,i+nat*(t-1))) )
                             ! 3D SISF
                             r2t_sisf_TMP(j)=r2t_sisf_TMP(j)+ cos(dot_product( ktd(:) , (pos(:,i)-pos0(:,i+nat*(t-1))) ))
                             ! Norm: - basically all the atoms in the 3D case. If not, something is fishy...
                             normz(j)=normz(j)+1
                          enddo
                       ! end 3D
                       endif
                    endif
                 enddo
                 ! 2D
                 ! Normalize these specific contributions and add it to the global arrays
                 if (td.ne.'yes') then
                    do j=1,n_zmesh
                       ! For some z-slices, it is possible we don't have atoms at
                       ! all. In that case, set r2t and r2t_sisf = 0 (instead of -nan...)
                       if (normz(j).eq.0) then
                          r2t(delt,j)=0.0
                          r2t_sisf(delt,j)=0.0
                       else      
                          r2t(delt,j)=r2t(delt,j)+r2t_TMP(j)/real(normz(j))
                          r2t_sisf(delt,j)=r2t_sisf(delt,j)+r2t_sisf_TMP(j)/real(normz(j))
                       endif
                       ! Square this guy, needed for chi4...
                       chi4(delt,j)=chi4(delt,j)+(r2t_sisf_TMP(j))**2.0 ! NO norm by N!
                       ! avoid norm N also for sisf -> needed for chi4
                       r2t_sisf_chi4(delt,j)=r2t_sisf_chi4(delt,j)+r2t_sisf_TMP(j) ! avoid norm by N!
                    enddo
                 ! end 2D
                 else
                 ! 3D
                    do j=1,n_zmesh
                       r2t(delt,j)=r2t(delt,j)+r2t_TMP(j)/real(normz(j))
                       r2t_sisf(delt,j)=r2t_sisf(delt,j)+r2t_sisf_TMP(j)/real(normz(j))
                       ! Square this guy, needed for chi4...
                       chi4(delt,j)=chi4(delt,j)+(r2t_sisf_TMP(j)/real(normz(j)))**2.0 ! normalized by N - this is 3D
                      ! ! avoid norm N also for sisf -> needed for chi4 in 2D
                      ! r2t_sisf_chi4(delt,j)=r2t_sisf_chi4(delt,j)+r2t_sisf_TMP(j) ! not norm by N! 
                    enddo
                 ! end 3D
                 endif
              endif
           enddo
       endif
       isamp=isamp+1
       if (isamp.gt.nframes) goto 1000
    endif
    STAT=read_xtc(xd,NATOMS,STEP,time,box_trans,pos,prec)
    counter=counter+1
enddo

STAT = xdrfile_close(xd)
deallocate(pos)

! check that we have specified the correct number of frames
1000 if (isamp-1.ne.nframes) then
        write(*,*) "Wrong number of frames! You said: ", nframes, " but I've read: ", isamp-1
     else
        write(*,*) isamp-1, " frames have been read"
        write(*,*) "We have selected ", ntel-1, " frames every ", nsamp, " frames"
        write(*,*) "We have selected ", t0-1, " time origins"
        write(*,*) "We have divided the system into ", n_zmesh, " slices, each one ", dz, " nm thick"
        if (td.eq.'yes') then
           write(*,*) "Note that this is a 3D calculation!!"
        endif
     endif

open(unit=102, file='hin_dynamics.out.tplots', status='unknown')
do j=1,n_zmesh
   write(102,*)
   write(102,*)
   do i=1,tmax-1
      time=dtime*(i-1)
      ! 2D
      if (td.ne.'yes') then
         ! We have already normalize each contribution from each time window by the - variable - number of atoms,
         ! all that is left is the normalization wrt the number of time windows...
         write(102,'(5f20.10)')time, r2t(i,j)/real(ntime(i)), r2t_sisf(i,j)/real(ntime(i)), &
              (chi4(i,j)/real(ntime(i)))-((r2t_sisf_chi4(i,j)/real(ntime(i)))**2.0), bound(j,1)+(dz/2.0) ! NO norm by N!
      ! end 2D
      else
      ! 3D  
         write(102,'(4f20.10)')time, r2t(i,j)/real(ntime(i)), r2t_sisf(i,j)/real(ntime(i)), &
              ((chi4(i,j)/real(ntime(i)))-((r2t_sisf(i,j)/real(ntime(i)))**2.0))*normz(1)
      endif
      ! 3D
   enddo
enddo

! Get the relaxation time, basically cutting the SISF at the value closest to
! TAUT(ttau) 
do j=1,n_zmesh
   ttau_min=1.0d30
   do i=1,tmax-1
      time=dtime*(i-1)
      tmin=abs(ttau-(r2t_sisf(i,j)/real(ntime(i)))) 
      if (tmin.lt.ttau_min) then
         ttau_min=tmin
         ttau_min_TIME=time      
      endif
   enddo
   ! Relaxation time is in ps
   relt(j)=ttau_min_TIME
enddo

! Get the diffusion coefficient as a function z
! We assume that the MSD is linear from (tmax-1)-((tmax-1)/linear)
tl_ini=(tmax-1)-int(real(tmax-1)/linear)
do j=1,n_zmesh
   n_linear=0
   do i=tl_ini,tmax-1
      n_linear=n_linear+1
      time=dtime*(i-1)
      sumx(j)=sumx(j)+time
      sumxsq(j)=sumxsq(j)+time**2.0 
      sumy(j)=sumy(j)+r2t(i,j)/real(ntime(i))
      sumxy(j)=sumxy(j)+time*r2t(i,j)/real(ntime(i))
   enddo
enddo
open(unit=103, file='hin_dynamics.out.dc_tau', status='unknown')
do j=1,n_zmesh
   ! Diffusion coefficient is in cm^2/s
   dc(j)=dcoeff*(((n_linear*sumxy(j))-(sumx(j)*sumy(j)))/((n_linear*sumxsq(j))-((sumx(j))**2.0)))
   write(103,'(3f20.10)') bound(j,1)+(dz/2.0), dc(j), relt(j)
enddo

end program hin_dynamics

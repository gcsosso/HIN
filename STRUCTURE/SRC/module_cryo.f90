module MOD_cryo

contains

subroutine cryo_alloc(pos,nat,sym,ns,n_ws,list_ws,o_ns,cart,icell,list_nw,n_nw,nr,dr,half_dr,rad,gr_norm,o_solv,n_solv)

implicit none

! Arguments
real :: icell(cart*cart), dr, half_dr
real, allocatable :: pos(:,:), rad(:), gr_norm(:,:)
character*4, allocatable :: sym(:)
integer :: nat, ns, cart, o_ns, n_nw, nr, n_solv
integer, allocatable :: n_ws(:), list_ws(:,:), list_nw(:), o_solv(:)

! Local
integer :: i, j, pairs, counter, ir
real :: l_box

do i=1,ns
   if (sym(list_ws(i,1)).eq.'OW') then
      o_ns=i ! o_ns = index of OW in list_ws
   endif
   !do j=1,n_ws(i)
   !   write(*,*) sym(list_ws(i,j))
   !enddo
enddo

n_nw=nat-n_ws(o_ns)*4 ! n_nw = number of non-water species. If using TIP4P water model i.e. 4 particles/molecule !! Hard-coded
allocate(list_nw(n_nw))

pairs=((n_ws(o_ns)-1)**2.0+(n_ws(o_ns)-1))/2.0 ! Number of O-O pairs

! Get indexes of all non-water species (list_nw)
counter=1
do i=1,nat
  if ((sym(i).eq.'OW').or.(sym(i).eq.'HW1').or.(sym(i).eq.'HW2').or.(sym(i).eq.'MW')) then ! Not sure if this is the best way to do this... could also search against a separate vector containing water species names
    continue
  else
    list_nw(counter)=i
    counter=counter+1
  endif
enddo

n_nw=size(list_nw) ! n_nw = number of non-water species

! Pair correlation functions: build mesh
l_box=icell(1)
allocate(rad(nr))
dr=l_box/(real(2.0*nr))
half_dr=dr/2.0

do ir=1,nr
  rad(ir)=(real(ir-1)+0.5)*dr
enddo

allocate(gr_norm(n_nw,nr),o_solv(n_ws(o_ns)))
gr_norm(:,:)=0.0d0
o_solv(:)=0
n_solv=0

end subroutine cryo_alloc


subroutine cryo(pos,sym,ns,n_ws,list_ws,o_ns,cart,icell,list_nw,n_nw,c_rcut,nr,dr,half_dr,rad,gr_norm,fact)

implicit none

! Arguments
real :: icell(cart*cart), c_rcut, dr, half_dr, fact
real, allocatable :: pos(:,:), rad(:), gr_norm(:,:)
character*4, allocatable :: sym(:)
integer :: ns, o_ns, cart, n_nw, nr
integer, allocatable :: n_ws(:), list_ws(:,:), list_nw(:)

! Local
integer :: i, j, counter, i_spc, j_spc, ir
integer, allocatable :: gr(:,:)
real :: i_pos(3), j_pos(3), xdf, ydf, zdf, r_ij, vol
real :: rho, r2, num_i, num_j, volume, gr_tmp
real(8), parameter :: pi=4.0*datan(1.d0), pi4=4.0*pi

allocate(gr(n_nw,nr))
gr(:,:)=0.0d0
gr_tmp=0.0d0

! O-O Pair correlation functions: count atoms, assign to bins
! do i=1,n_ws(o_ns)
!   i_spc=list_ws(o_ns,i)
!   i_pos(1)=pos(1,i_spc) ; i_pos(2)=pos(2,i_spc) ; i_pos(3)=pos(3,i_spc)
!
!   do j=1,n_ws(o_ns)
!     j_spc=list_ws(o_ns,j)
!     if (i_spc.eq.j_spc) cycle
!     j_pos(1)=pos(1,j_spc) ; j_pos(2)=pos(2,j_spc) ; j_pos(3)=pos(3,j_spc)
!
!     xdf=i_pos(1)-j_pos(1)
!     ydf=i_pos(2)-j_pos(2)
!     zdf=i_pos(3)-j_pos(3)
!
!     call images(cart,0,1,1,icell,xdf,ydf,zdf)
!     r_ij=sqrt(xdf**2.0+ydf**2.0+zdf**2.0)
!
!
!     do ir=1,nr
!       if ((r_ij.gt.rad(ir)-half_dr).and.(r_ij.le.rad(ir)+half_dr)) then
!          gr(ir)=gr(ir)+1
!       endif
!     enddo
!   enddo
! enddo

! Mol-O Pair correlation functions: count atoms, assign to bins
do i=1,n_nw
  i_spc=list_nw(i)
  i_pos(1)=pos(1,i_spc) ; i_pos(2)=pos(2,i_spc) ; i_pos(3)=pos(3,i_spc)

  do j=1,n_ws(o_ns)
    j_spc=list_ws(o_ns,j)
    j_pos(1)=pos(1,j_spc) ; j_pos(2)=pos(2,j_spc) ; j_pos(3)=pos(3,j_spc)

    xdf=i_pos(1)-j_pos(1)
    ydf=i_pos(2)-j_pos(2)
    zdf=i_pos(3)-j_pos(3)

    call images(cart,0,1,1,icell,xdf,ydf,zdf)
    r_ij=sqrt(xdf**2.0d0+ydf**2.0d0+zdf**2.0d0)

    do ir=1,nr
      if ((r_ij.gt.rad(ir)-half_dr).and.(r_ij.le.rad(ir)+half_dr)) then
        gr(i,ir)=gr(i,ir)+1
      endif
    enddo
  enddo
enddo

! O-O PCF normalisation
! num_i=dble(n_ws(o_ns))
! num_j=num_i
! volume=icell(1)**3.0d0

! Mol-O PCF normalisation
!num_i=dble(n_nw) ! use if calculating a PCF for whole molecule
num_i=1.0d0 ! use if calculating a PCF for individual atoms
num_j=dble(n_ws(o_ns))
volume=icell(1)**3.0d0

fact=pi4*dr*(num_j/volume)
do i=1,n_nw
  do ir=1,nr
    r2=(rad(ir))**2.0d0
    gr_tmp=0.0d0 ! Reset
    gr_tmp=gr(i,ir)/(fact*r2*num_i)
    gr_norm(i,ir)=gr_norm(i,ir)+gr_tmp
    !write(88,*) rad(ir), gr_norm(ir)
  enddo
enddo

end subroutine cryo


subroutine cryo_workup(fframe,lframe,n_nw,rad,dr,gr_norm,gr_average,smgr_average,cn_running,rmin)

! Arguments
real :: dr
real, allocatable :: rad(:), gr_norm(:,:), gr_average(:,:), smgr_average(:,:), cn_running(:,:), rmin(:)
integer :: fframe, lframe, n_nw

! Local
real :: density
integer :: i, j, k, n_bins, width
logical :: min

allocate(gr_average(n_nw,n_bins), smgr_average(n_nw,n_bins), cn_running(n_nw,n_bins), rmin(n_nw))
n_bins=size(rad)

write(*,*) "DEBUG", n_nw, n_bins

! gr (averaged over the n. of frames)
gr_average(:,:)=0.0d0
do i=1,n_nw
  do j=1,n_bins
    gr_average(i,j)=gr_norm(i,j)/dble(lframe-fframe+1)
  enddo
enddo

! Running coordination number
cn_running(:,:)=0.0d0
density=34.34375 !! hardcoded...
do i=1,n_nw
  do j=1,n_bins
   cn_running(i,j)=cn_running(i,j)+(gr_average(i,j)*rad(j)*rad(j)*dr) ! integration
 enddo
enddo

cn_running(:,:)=cn_running(:,:)*4.0d0*2.D0*DASIN(1.D0)*density ! normalisation

! smooth
smgr_average(:,:)=0.0d0
do i=1,n_nw
  do j=1,n_bins
    if (j.le.2) then
      smgr_average(i,j)=gr_average(i,3)
    elseif (j.ge.n_bins-2) then
      smgr_average(i,j)=gr_average(i,n_bins)
    else
      smgr_average(i,j)=(gr_average(i,j-2)+2.0d0*gr_average(i,j-1)+3.0d0*gr_average(i,j)+2.0d0*gr_average(i,j+1)+gr_average(i,j+2))/9.0d0 !! Hard coded - use should be able to enter smoothing coarseness
    endif
  enddo
enddo

! Find first minimum with search width parameter
rmin(:)=0.0d0
width=2 !! Hard coded
min=.false.
do i=1,n_nw
  do j=1,n_bins
    if (j.le.width .or. j.ge.n_bins-width) then
      cycle ! skip - otherwise will go out of bounds at next conditional statement
    else
      do k=1, width
        if (smgr_average(i,j-k).gt.smgr_average(i,j) .and. smgr_average(i,j+k).gt.smgr_average(i,j)) then
          min=.true.
        else
          min=.false.
          exit ! rad(i) not the minimum - exit the inner loop
        endif
      enddo
      if (min.eqv..true.) then
        rmin(i)=rad(j)
        min=.false.
        exit ! minimum found - exit middle loop
      endif
    endif
  enddo
enddo

end subroutine cryo_workup

subroutine hydration(pos,n_ws,list_ws,o_ns,list_nw,n_nw,o_solv,n_solv,rmin,icell)

! Local
real :: i_pos(3), j_pos(3), xdf, ydf, zdf, r_ij
integer :: i, j, k, i_spc, j_spc
integer, parameter :: cart=3

! Arguments
real, allocatable :: pos(:,:), rmin(:)
integer :: o_ns, n_nw
integer, allocatable :: n_ws(:), list_ws(:,:), list_nw(:), o_solv(:)
real :: icell(cart*cart)

!do i=1,n_nw
!   write(*,*) "DEBUG", rmin(i)
!enddo

do i=1,n_nw
  i_spc=list_nw(i)
  i_pos(1)=pos(1,i_spc) ; i_pos(2)=pos(2,i_spc) ; i_pos(3)=pos(3,i_spc)

  do j=1,n_ws(o_ns)
    j_spc=list_ws(o_ns,j)
    j_pos(1)=pos(1,j_spc) ; j_pos(2)=pos(2,j_spc) ; j_pos(3)=pos(3,j_spc)

    xdf=i_pos(1)-j_pos(1)
    ydf=i_pos(2)-j_pos(2)
    zdf=i_pos(3)-j_pos(3)

    !write(*,*) "DEBUG_1", i, cart, icell(:), xdf, ydf, zdf 
    call images(cart,0,1,1,icell,xdf,ydf,zdf)
    r_ij=sqrt(xdf**2.0d0+ydf**2.0d0+zdf**2.0d0)
    !write(*,*) "DEBUG_2",i

    if (r_ij.lt.rmin(i)) then
        o_solv(j)=1
    endif
  enddo
enddo

n_solv=n_solv+sum(o_solv)

end subroutine hydration

end module MOD_cryo

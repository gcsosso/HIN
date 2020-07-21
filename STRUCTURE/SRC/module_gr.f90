module MOD_gr

contains

subroutine gr_alloc(nat,sym,ns,n_ws,list_ws,o_ns,cart,icell,list_nw,n_nw,n_ow,gr_bins,dr,half_dr,rad,gr_mol_norm,gr_atm_norm,o_dist)

implicit none

! Arguments
real :: icell(cart*cart), dr, half_dr
real, allocatable :: rad(:), gr_mol_norm(:), gr_atm_norm(:,:), o_dist(:)
character*4, allocatable :: sym(:)
integer :: nat, ns, cart, o_ns, n_nw, n_ow, gr_bins
integer, allocatable :: n_ws(:), list_ws(:,:), list_nw(:)

! Local
integer :: i, n_oo, counter
real :: l_box

! Get list of OW indices from list_ws !! Will "OW" work for all water models??
do i=1,ns
   if (sym(list_ws(i,1)).eq.'OW') then
      o_ns=i
   endif
enddo

! Get number of non-water species to allocate list_nw. !! Hard-coded for TIP4P (*4 for 4 particles in model)
n_nw=nat-n_ws(o_ns)*4
allocate(list_nw(n_nw))

! Number of Ow species and O-O pairs
n_ow=n_ws(o_ns)
n_oo=((n_ow-1)**2.0+(n_ow-1))/2.0

! Get indexes of all non-water species (list_nw) !! Also hard-coded for TIP4P particle names
counter=1
do i=1,nat
  if ((sym(i).eq.'OW').or.(sym(i).eq.'HW1').or.(sym(i).eq.'HW2').or.(sym(i).eq.'MW')) then
    continue
  else
    list_nw(counter)=i
    counter=counter+1
  endif
enddo

! Build mesh
l_box=icell(1)
dr=l_box/(real(2.0*gr_bins))
half_dr=dr/2.0

allocate(rad(gr_bins))
do i=1,gr_bins
  rad(i)=(real(i-1)+0.5)*dr
enddo

allocate(gr_mol_norm(gr_bins),gr_atm_norm(n_nw,gr_bins),o_dist(n_ow))
gr_mol_norm(:)=0.0d0
gr_atm_norm(:,:)=0.0d0

end subroutine gr_alloc


subroutine gr(pos,list_ws,o_ns,cart,icell,list_nw,n_nw,n_ow,gr_bins,dr,half_dr,rad,gr_mol_norm,gr_atm_norm,fact,o_dist)

implicit none

! Arguments
real :: icell(cart*cart), dr, half_dr, fact
real, allocatable :: pos(:,:), rad(:), gr_mol_norm(:), gr_atm_norm(:,:), o_dist(:)
integer :: o_ns, cart, n_nw, n_ow, gr_bins
integer, allocatable :: list_ws(:,:), list_nw(:)

! Local
integer :: i, j, i_spc, j_spc, ir
integer, allocatable :: gr_mol(:), gr_atm(:,:)
real :: i_pos(3), j_pos(3), xdf, ydf, zdf, r_ij
real :: r2, num_i, num_j, volume, gr_tmp
real(8), parameter :: pi=4.0*datan(1.d0), pi4=4.0*pi

allocate(gr_mol(gr_bins),gr_atm(n_nw,gr_bins))
gr_mol(:)=0
gr_atm(:,:)=0
gr_tmp=0.0d0

! O-O PCF [I]
! do i=1,n_ow
!   i_spc=list_ws(o_ns,i)
!   i_pos(1)=pos(1,i_spc) ; i_pos(2)=pos(2,i_spc) ; i_pos(3)=pos(3,i_spc)
!
!   do j=1,n_ow
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
!     do ir=1,gr_bins
!       if ((r_ij.gt.rad(ir)-half_dr).and.(r_ij.le.rad(ir)+half_dr)) then
!          gr_mol(ir)=gr_mol(ir)+1
!       endif
!     enddo
!   enddo
! enddo

! M-O PCF: individual atoms [II]
do i=1,n_nw
  i_spc=list_nw(i)
  i_pos(1)=pos(1,i_spc) ; i_pos(2)=pos(2,i_spc) ; i_pos(3)=pos(3,i_spc)

  do j=1,n_ow
    j_spc=list_ws(o_ns,j)
    j_pos(1)=pos(1,j_spc) ; j_pos(2)=pos(2,j_spc) ; j_pos(3)=pos(3,j_spc)

    xdf=i_pos(1)-j_pos(1)
    ydf=i_pos(2)-j_pos(2)
    zdf=i_pos(3)-j_pos(3)

    call images(cart,0,1,1,icell,xdf,ydf,zdf)
    r_ij=sqrt(xdf**2.0d0+ydf**2.0d0+zdf**2.0d0)

    do ir=1,gr_bins
      if ((r_ij.gt.rad(ir)-half_dr).and.(r_ij.le.rad(ir)+half_dr)) then
        gr_atm(i,ir)=gr_atm(i,ir)+1
      endif
    enddo
  enddo
enddo

! M-O PCF: whole molecule, then averaged over atoms (to account for multiple counting) [III]
! do i=1,n_nw
!   i_spc=list_nw(i)
!   i_pos(1)=pos(1,i_spc) ; i_pos(2)=pos(2,i_spc) ; i_pos(3)=pos(3,i_spc)
!
!   do j=1,n_ow
!     j_spc=list_ws(o_ns,j)
!     j_pos(1)=pos(1,j_spc) ; j_pos(2)=pos(2,j_spc) ; j_pos(3)=pos(3,j_spc)
!
!     xdf=i_pos(1)-j_pos(1)
!     ydf=i_pos(2)-j_pos(2)
!     zdf=i_pos(3)-j_pos(3)
!
!     call images(cart,0,1,1,icell,xdf,ydf,zdf)
!     r_ij=sqrt(xdf**2.0d0+ydf**2.0d0+zdf**2.0d0)
!
!     do ir=1,gr_bins
!       if ((r_ij.gt.rad(ir)-half_dr).and.(r_ij.le.rad(ir)+half_dr)) then
!         gr_mol(ir)=gr_mol(ir)+1
!       endif
!     enddo
!   enddo
! enddo

! M-O PCF: whole molecule, treated as single particle (using smallest M-Ow distances) [IV]
! do i=1,n_nw
!   i_spc=list_nw(i)
!   i_pos(1)=pos(1,i_spc) ; i_pos(2)=pos(2,i_spc) ; i_pos(3)=pos(3,i_spc)
!
!   do j=1,n_ow
!     j_spc=list_ws(o_ns,j)
!     j_pos(1)=pos(1,j_spc) ; j_pos(2)=pos(2,j_spc) ; j_pos(3)=pos(3,j_spc)
!
!     xdf=i_pos(1)-j_pos(1)
!     ydf=i_pos(2)-j_pos(2)
!     zdf=i_pos(3)-j_pos(3)
!
!     call images(cart,0,1,1,icell,xdf,ydf,zdf)
!     r_ij=sqrt(xdf**2.0d0+ydf**2.0d0+zdf**2.0d0)
!
!     if (i.eq.1) then ! First iteration of outer loop, assign regardless of value
!       o_dist(j)=r_ij
!     elseif (r_ij .lt. o_dist(j)) then ! If the water is closer to another atom, assign this new shorter distance
!       o_dist(j)=r_ij
!     endif
!   enddo
! enddo
! do i=1,n_ow
!   do j=1,gr_bins
!     if ((o_dist(i).gt.rad(j)-half_dr).and.(o_dist(i).le.rad(j)+half_dr)) then
!       gr_mol(j)=gr_mol(j)+1
!     endif
!   enddo
! enddo

! Normalisation
! num_i=dble(n_ow) ! Use for O-O PCF [I]
! num_i=dble(n_nw) ! Use for M-O PCF (averaged over atoms [III])
num_i=1.0d0 ! Use for M-O PCF (individual atoms [II] or as single particle [IV])
num_j=dble(n_ow)
volume=icell(1)**3.0d0
fact=pi4*dr*(num_j/volume)

! Use for M-O PCF (individual atoms [II])
do i=1,n_nw
  do ir=1,gr_bins
    r2=(rad(ir))**2.0d0
    gr_tmp=gr_atm(i,ir)/(fact*r2*num_i)
    gr_atm_norm(i,ir)=gr_atm_norm(i,ir)+gr_tmp
  enddo
enddo

! Use for O-O PCF [I] or M-O PCF (averaged [III] or minimum distance [IV])
! do ir=1,gr_bins
!   r2=(rad(ir))**2.0d0
!   gr_tmp=gr_mol(ir)/(fact*r2*num_i)
!   gr_mol_norm(ir)=gr_mol_norm(ir)+gr_tmp
! enddo

end subroutine gr

end module MOD_gr

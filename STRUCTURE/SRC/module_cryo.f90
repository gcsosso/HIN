module MOD_cryo

contains

subroutine cryo_alloc(pos,nat,sym,ns,n_ws,list_ws,o_dist,o_ns,list_nw,n_nw,o_rad_count,nr,dr,half_dr,rad)

implicit none

! Arguments
real :: o_rad_count, dr, half_dr
real, allocatable :: pos(:,:), o_dist(:), rad(:)
character*4, allocatable :: sym(:)
integer :: nat, ns, o_ns, n_nw, nr
integer, allocatable :: n_ws(:), list_ws(:,:), list_nw(:)

! Local
integer :: i, j, pairs, counter, ir
real :: l_box

o_rad_count=0

do i=1,ns
   if (sym(list_ws(i,1)).eq.'OW') then
      o_ns=i ! o_ns = index of OW in list_ws
   endif
   !do j=1,n_ws(i)
   !   write(*,*) sym(list_ws(i,j))
   !enddo
enddo

n_nw=nat-n_ws(o_ns)*4 ! n_nw = number of non-water species. If using TIP4P water model i.e. 4 particles/molecule
allocate(list_nw(n_nw))

pairs=((n_ws(o_ns)-1)**2.0+(n_ws(o_ns)-1))/2.0 ! Number of O-O pairs
allocate(o_dist(pairs))

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
l_box=4.0 ! TBE
allocate(rad(nr))
dr=l_box/(real(2.0*nr))
half_dr=dr/2.0

do ir=1,nr
  rad(ir)=(real(ir-1)+0.5)*dr
enddo

end subroutine cryo_alloc


subroutine cryo(pos,sym,ns,n_ws,list_ws,o_dist,o_ns,cart,icell,list_nw,n_nw,c_rcut,o_rad_count,nr,dr,half_dr,rad)

implicit none

! Arguments
real :: icell(cart*cart), c_rcut, o_rad_count, dr, half_dr
real, allocatable :: pos(:,:), o_dist(:), rad(:)
character*4, allocatable :: sym(:)
integer :: ns, o_ns, cart, n_nw, nr
integer, allocatable :: n_ws(:), list_ws(:,:), list_nw(:)

! Local
integer :: i, j, counter, i_spc, j_spc, ir
integer, allocatable :: o_rad(:)
real :: i_pos(3), j_pos(3), xdf, ydf, zdf, r_ij, vol, rho, r2, gr_norm
real, allocatable :: gr(:,:)
real(8), parameter :: pi=4.0*datan(1.d0)

! Calculate distances between each pair of atoms in o_ns
! counter=1
! do i=1,n_ws(o_ns)-1
!   i_pos(1)=pos(1,list_ws(o_ns,i))
!   i_pos(2)=pos(2,list_ws(o_ns,i))
!   i_pos(3)=pos(3,list_ws(o_ns,i))
!
!   do j=i+1,n_ws(o_ns)
!       j_pos(1)=pos(1,list_ws(o_ns,j))
!       j_pos(2)=pos(2,list_ws(o_ns,j))
!       j_pos(3)=pos(3,list_ws(o_ns,j))
!
!       xdf=i_pos(1)-j_pos(1)
!       ydf=i_pos(2)-j_pos(2)
!       zdf=i_pos(3)-j_pos(3)
!
!       call images(cart,0,1,1,icell,xdf,ydf,zdf) ! Images subroutine returns xdf, ydf and zdf updated to account for PBCs
!       ij_dist=sqrt(xdf**2.0+ydf**2.0+zdf**2.0)
!       o_dist(counter)=ij_dist
!       counter=counter+1
!    enddo
! enddo

allocate(o_rad(n_ws(o_ns)), gr(n_nw,nr))
o_rad=0 ! Vector of 0s corresponding to O water atoms - to be 'coloured in' when atoms are counted

! Pair correlation functions: count atoms, assign to bins
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
    r_ij=sqrt(xdf**2.0+ydf**2.0+zdf**2.0)

    !if (r_ij.lt.c_rcut) then
      !o_rad(j)=1 ! 'colour in' the corresponding cell

    do ir=1,nr
      if ((r_ij.gt.rad(ir)-half_dr).and.(r_ij.le.rad(ir)+half_dr)) then
        gr(i,ir)=gr(i,ir)+1
      endif
    enddo
  enddo
enddo

! Pair correlation functions: normalise
vol=4.0**3 ! TBE
rho=real(n_ws(o_ns)+n_nw)/vol ! Should this be density of water oxygens or water oxygens and molecule of interest?

do ir=1,nr
  r2=rad(ir)**2

  do i=1,n_nw
    !gr(i,ir)=gr(i,ir)/(4.0*pi*r2*dr*rho*(n_nw+n_ws(o_ns)))
    gr_norm=gr(i,ir)/n_ws(o_ns) ! Number of particles (water Os) considered
    gr_norm=gr_norm/(4.0*pi*r2*dr) ! Volume of spherical shell
    gr_norm=gr_norm/rho ! Particle density
    gr(i,ir)=gr_norm

  enddo
enddo

! o_rad_count=o_rad_count+sum(o_rad)
! write(*,*) o_rad_count

! ! Pair correlation functions: output
! open(unit=19,file="cryo_out.dat",status="new")
! do i=1,n_nw
!   write(19,*) sym(list_nw(i))
!   do ir=1,nr
!     write(19,*) rad(ir), gr(i,ir)
!   enddo
!   write(19,*)
! enddo

! Debug
write(*,*) sym(list_nw(2)), dr
do ir=1,nr
    write(*,*) rad(ir), gr(2,ir)
enddo

end subroutine cryo

end module MOD_cryo

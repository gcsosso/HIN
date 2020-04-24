module MOD_cryo

contains

subroutine cryo_alloc(pos,nat,sym,ns,n_ws,list_ws,o_dist,o_ns,list_nw,n_nw)

implicit none

! Arguments
real, allocatable :: pos(:,:), o_dist(:)
character*4, allocatable :: sym(:)
integer :: nat, ns, o_ns, n_nw
integer, allocatable :: n_ws(:), list_ws(:,:), list_nw(:)

! Local
integer :: i, j, pairs, counter

do i=1,ns
   if (sym(list_ws(i,1)).eq.'OW') then
      o_ns=i ! o_ns = index of OW in list_ws
   endif
   !do j=1,n_ws(i)
   !   write(*,*) sym(list_ws(i,j))
   !enddo
enddo

n_nw=nat-n_ws(o_ns)*4.0 ! n_nw = number of non-water species. If using TIP4P water model i.e. 4 particles/molecule
allocate(list_nw(n_nw))

pairs=((n_ws(o_ns)-1)**2.0+(n_ws(o_ns)-1))/2.0 ! Number of O-O pairs
allocate(o_dist(pairs))

! Get indexes of all non-water species (list_nw)
counter=1
do i=1,nat
  if ( (sym(i).eq.'OW') .OR. (sym(i).eq.'HW1') .OR. (sym(i).eq.'HW2') .OR. (sym(i).eq.'MW') ) then ! Not sure if this is the best way to do this... could also search against a separate vector containing water species names
    continue
  else
    list_nw(counter)=i
    counter=counter+1
  endif
enddo

n_nw=size(list_nw) ! n_nw = number of non-water species

end subroutine cryo_alloc

subroutine cryo(pos,sym,ns,n_ws,list_ws,o_dist,o_ns,cart,icell,list_nw,n_nw,c_rcut)

implicit none

! Arguments
real :: icell(cart*cart), c_rcut
real, allocatable :: pos(:,:), o_dist(:)
character*4, allocatable :: sym(:)
integer :: ns, o_ns, cart, n_nw
integer, allocatable :: n_ws(:), list_ws(:,:), list_nw(:)

! Local
integer :: i, j, counter
integer, allocatable :: o_rad(:)
real :: i_pos(3), j_pos(3), xdf, ydf, zdf, ij_dist

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

allocate(o_rad(n_ws(o_ns)))
o_rad = 0

! Calculate distances between each non-water atom (i) and oxygen (j)
do i=1,n_nw
  i_pos(1)=pos(1,list_nw(i))
  i_pos(2)=pos(2,list_nw(i))
  i_pos(3)=pos(3,list_nw(i))

  do j=1,n_ws(o_ns)
    j_pos(1)=pos(1,list_ws(o_ns,j))
    j_pos(2)=pos(2,list_ws(o_ns,j))
    j_pos(3)=pos(3,list_ws(o_ns,j))

    xdf=i_pos(1)-j_pos(1)
    ydf=i_pos(2)-j_pos(2)
    zdf=i_pos(3)-j_pos(3)

    call images(cart,0,1,1,icell,xdf,ydf,zdf)
    ij_dist=sqrt(xdf**2.0+ydf**2.0+zdf**2.0)

    if (ij_dist.lt.c_rcut) then
      o_rad(j)=1 ! 'colour in' the corresponding cell
    endif
  enddo
enddo

write(*,*) sum(o_rad)

! Next steps:
! - The number of O atoms within the cutoff (sum(o_rad)) needs to be stored for each frame so that an average can be computed for the entire trajectory
! - Different cutoff radii should be included for different species. These should be specified in the input file
!     - We can do this by including a input line for all species found in molecule of interest and values for each

end subroutine cryo

end module MOD_cryo

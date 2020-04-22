module MOD_cryo

contains

! Usually, in the HIN code, each module contains two subroutines:
! 1. An "alloc" subroutine.
! 2. A subroutine that actually computes stuff

subroutine cryo_alloc(pos,sym,ns,n_ws,list_ws,o_dist,o_ns)

implicit none

! Arguments
real, allocatable :: pos(:,:), o_dist(:)
character*4, allocatable :: sym(:)
integer :: ns, o_ns
integer, allocatable :: n_ws(:), list_ws(:,:)

! Local
integer :: i, j, pairs

do i=1,ns
   if (sym(list_ws(i,1)).eq.'OW') then
      !write(*,*) "ns=i -> these are the oxygens", i
      o_ns=i
   endif
   !do j=1,n_ws(i)
   !   write(*,*) sym(list_ws(i,j))
   !enddo
enddo

pairs=((n_ws(o_ns)-1)**2+(n_ws(o_ns)-1))/2
!write(*,*) n_ws(o_ns), pairs

allocate(o_dist(pairs))

end subroutine cryo_alloc


subroutine cryo(pos,sym,ns,n_ws,list_ws,o_dist,o_ns,cart,icell)

implicit none

! Arguments
real :: icell(cart*cart)
real, allocatable :: pos(:,:), o_dist(:)
character*4, allocatable :: sym(:)
integer :: ns, o_ns, cart
integer, allocatable :: n_ws(:), list_ws(:,:)

! Local variables
integer :: i, j, k, counter
real :: x1, y1, z1, x2, y2, z2, xdf, ydf, zdf, ij_dist

!do i=1,ns
!   write(*,*) ns, n_ws(i)
!enddo

! Calculate distances between each pair of atoms in o_ns
! cart = n. of cartesian directions (3)
! 0, 1, 1 - do not touch them!
! icell: 1D vector containing 9 elements - the elements of the cell matrix. -> pass icell into the subroutine
! xdf, ydf, zdf = x2-x1, y2-y1, z2-z1
! The images subroutine returns the SAME xdf, ydf and zdf (they get overwritten) BUT this time modified so as to
! take into account the PBCs

counter=1
do i=1,n_ws(o_ns)-1

  x1=pos(1,list_ws(o_ns,i))
  y1=pos(2,list_ws(o_ns,i))
  z1=pos(3,list_ws(o_ns,i))

  do j=i+1,n_ws(o_ns)

      x2=pos(1,list_ws(o_ns,j))
      y2=pos(2,list_ws(o_ns,j))
      z2=pos(3,list_ws(o_ns,j))

      xdf=x1-x2 ; ydf=y1-y2 ; zdf=z1-z2

      call images(cart,0,1,1,icell,xdf,ydf,zdf) ! Images subroutine returns xdf, ydf and zdf updated to account for PBCs

      ij_dist=sqrt(xdf**2+ydf**2+zdf**2)
      o_dist(counter)=ij_dist

      counter=counter+1

   enddo
enddo

! Lets test this out
do k=1,3
  write(*,*) o_dist(k)
enddo

end subroutine cryo

end module MOD_cryo

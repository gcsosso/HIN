module MOD_cryo

contains

! Usually, in the HIN code, each module contains two subroutines:
! 1. An "alloc" subroutine.
! 2. A subroutine that actually computes stuff

subroutine cryo_alloc(pos,o_dist,sym,ns,n_ws,list_ws,ox_ns)

implicit none

! Arguments
real, allocatable :: pos(:,:), o_dist(:)
character*4, allocatable :: sym(:)
integer :: ns, ox_ns
integer, allocatable :: n_ws(:), list_ws(:,:)

! Local
integer :: i, j, pairs

do i=1,ns
   if (sym(list_ws(i,1)).eq.'OW') then
      !write(*,*) "ns=i -> these are the oxygens", i      
      ox_ns=i
   endif
   !do j=1,n_ws(i)
   !   write(*,*) sym(list_ws(i,j))
   !enddo
enddo

pairs=( ( n_ws(ox_ns)-1 )**2 + ( n_ws(ox_ns)-1) )/2
!write(*,*) n_ws(ox_ns), pairs

allocate(o_dist(pairs))

end subroutine cryo_alloc

! 

subroutine cryo(pos,sym,ns,n_ws,list_ws,o_dist,ox_ns)

implicit none

! Arguments
real, allocatable :: pos(:,:), o_dist(:)
character*4, allocatable :: sym(:)
integer :: ns, ox_ns
integer, allocatable :: n_ws(:), list_ws(:,:)

! Local variables
integer :: i, j, counter
real :: x1, y1, z1, x2, y2, z2

!do i=1,ns
!   write(*,*) ns, n_ws(i)
!enddo

!! Calculate distances between each pair of atoms in matrix
counter=1
do i=1,n_ws(ox_ns)-1
   do j=i+1,n_ws(ox_ns)
      x1=pos(1,list_ws(ox_ns,i)) 
      y1=pos(2,list_ws(ox_ns,i))
      z1=pos(3,list_ws(ox_ns,i))
      x2=pos(1,list_ws(ox_ns,j))
      y2=pos(2,list_ws(ox_ns,j))
      z2=pos(3,list_ws(ox_ns,j))

      ! We have the components of the positions of the oxygen atoms we are interested in.
      ! We need to find their distances - taking into account PBCs
      ! To this end we can use the "images" subroutine
      ! See e.g. module_rings

      ! call images(cart,0,1,1,icell,xdf,ydf,zdf)
      ! cart = n. of cartesian directions (3)
      ! 0, 1, 1 - do not touch them!
      ! icell: 1D vector containing 9 elements - the elements of the cell matrix. -> pass icell into the subroutine
      ! xdf, ydf, zdf = x2-x1, y2-y1, z2-z1 
      ! The images subroutine returns the SAME xdf, ydf and zdf (they get overwritten) BUT this time modified so as to
      ! take into account the PBCs 

      !±dist=sqrt(((x1-x2)**2)+((y1-y2)**2)+((z1-z2)**2))
      !o_dist(counter)=dist
      counter=counter+1
   enddo
enddo

end subroutine cryo


end module MOD_cryo

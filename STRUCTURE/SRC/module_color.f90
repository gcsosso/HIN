module MOD_color

contains

subroutine icy_water(lframe,idx_cls)
  implicit none
  integer :: i, lframe
  integer, allocatable :: idx_cls(:)
! read in icy OW for SOLVATION module
  open(unit=100, file='./idx.dat', status='old')
  allocate(idx_cls(lframe))
  do i=1,lframe
        read(100,*) idx_cls(i)
        !write(*,*) idx_cls(i)
  enddo
end subroutine icy_water

subroutine color_OW(idx_cls,counter,nat,coloring)
   integer :: j,nat,counter
   integer, allocatable :: icy(:),coloring(:),idx_cls(:)
   allocate(icy(idx_cls(counter+1)))
   !allocate(coloring(nat))
   coloring(:) = 0
   read(100,*) icy(:)

   do j=1,idx_cls(counter+1)
      coloring(icy(j)+1)=1
   enddo

   deallocate(icy)
end subroutine color_OW

end module MOD_color

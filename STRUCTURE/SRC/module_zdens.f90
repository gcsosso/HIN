module MOD_zdens

contains

subroutine zdens_alloc(nz,zmax,zmin,dz,dens,zmesh,ns)

implicit none

! Local
integer :: k

! Arguments
integer :: nz, ns
real :: zmax, zmin, dz
real, allocatable :: dens(:,:), zmesh(:)

! Allocate the z-dens array...
nz=int((zmax-zmin)/dz)
if (nz.lt.1) then
   write(99,*) "Your mesh along z (z-dens section) is ill-defined. Modify zmin, zmax and/or dz"
   stop
endif
allocate(dens(nz,ns),zmesh(nz))
dens(:,:)=0.0d0
do k=1,nz
   zmesh(k)=zmin+((k*dz)-(dz/2.0d0))
enddo

end subroutine zdens_alloc

subroutine zdens(ns,n_ws,nz,zmin,dz,pos,cart,list_ws,dens)

implicit none

! Local
integer :: i, j, k
real :: lb, ub

! Arguments
integer :: ns, nz, cart
integer, allocatable :: n_ws(:), list_ws(:,:)
real :: zmin, dz
real, allocatable :: pos(:,:), dens(:,:)

character*3 :: switch_zdens

! Number density profile along z...
do i=1,ns
   do j=1,n_ws(i)
      do k=1,nz
         lb=zmin+(k-1)*dz
         ub=lb+dz
         if (pos(cart,list_ws(i,j)).gt.lb.and.pos(cart,list_ws(i,j)).le.ub) then
            dens(k,i)=dens(k,i)+1.0d0
         endif
      enddo
   enddo
enddo

return

end subroutine zdens

end module MOD_zdens

module MOD_xyfes 

contains

subroutine xyfes_alloc(nxy,xymax,xymin,ddx,ddy,xydens,xmesh,ymesh,ns,icell,cart)

implicit none

! Local
integer :: k, l

! Arguments
integer :: nxy, ns, cart
real :: xymax, xymin, ddx, ddy, icell(cart*cart)
real, allocatable :: xydens(:,:,:), xmesh(:), ymesh(:)

! Allocate the xydens array...
ddx=(icell(1)/real(nxy))
ddy=(icell(5)/real(nxy))

allocate(xydens(nxy,nxy,ns),xmesh(nxy),ymesh(nxy))
xydens(:,:,:)=0.0d0

do k=1,nxy
   xmesh(k)=((k*ddx)-(ddx/2.0d0))
   ymesh(k)=((k*ddy)-(ddy/2.0d0))
enddo

end subroutine xyfes_alloc

subroutine xyfes(ns,n_ws,nxy,xymin,ddx,ddy,pos,cart,list_ws,xydens,xymax,icell)

implicit none

! Local
integer :: i, j, k, l, hww
real :: lbx, ubx, lby, uby, tmp_x, tmp_y

! Arguments
integer :: ns, nxy, cart
integer, allocatable :: n_ws(:), list_ws(:,:)
real :: xymin, ddx, ddy, xymax, icell(cart*cart)
real, allocatable :: pos(:,:), xydens(:,:,:)

! 2D FES in the xy plane...
do i=1,ns
   hww=0
   do j=1,n_ws(i)
      ! Select only those atoms within the selected slice...
      if (pos(cart,list_ws(i,j)).gt.xymin.and.pos(cart,list_ws(i,j)).lt.xymax) then
         hww=hww+1
         do k=1,nxy
            lbx=(k-1)*ddx
            ubx=lbx+ddx
            do l=1,nxy
               lby=(l-1)*ddy
               uby=lby+ddy  
               !write(*,*) j, pos(1,list_ws(i,j)), lbx, ubx, pos(2,list_ws(i,j)), lby, uby
               ! Orthorhombic PBCs...
               !if (list_ws(i,j).eq.126) write(*,*) "before", pos(1,list_ws(i,j)), pos(2,list_ws(i,j))
               tmp_x=pos(1,list_ws(i,j)) ; tmp_y=pos(2,list_ws(i,j))
               if (tmp_x.lt.0.0d0)             tmp_x=tmp_x-(icell(1)*floor(tmp_x/icell(1))) 
               if (tmp_x.ge.icell(1))          tmp_x=tmp_x+(icell(1)*floor(tmp_x/icell(1)))
               if (tmp_y.lt.0.0d0)             tmp_y=tmp_y-(icell(5)*floor(tmp_y/icell(5)))
               if (tmp_y.ge.icell(5))          tmp_y=tmp_y+(icell(5)*floor(tmp_y/icell(5)))
               !if (list_ws(i,j).eq.126) write(*,*) "after", tmp_x, tmp_y
               if (tmp_x.gt.lbx.and. &
                   tmp_x.le.ubx.and. &
                   tmp_y.gt.lby.and. &
                   tmp_y.le.uby) then
                  xydens(k,l,i)=xydens(k,l,i)+1.0d0
               endif 
            enddo
         enddo
      endif
   enddo
   !if (hww.eq.0) then
   !   write(*,*) "ACTHUNG - No atoms within the selected xymin and xymax boundaries!" 
   !else
   !   write(*,*) hww
   !endif
enddo

return

end subroutine xyfes

end module MOD_xyfes

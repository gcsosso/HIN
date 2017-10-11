module MOD_electro

contains

subroutine electro_alloc(e_zmin,e_zmax,e_dz,d_charge,e_zmesh, &
           e_nz,atq,qqq,nq,qqq_all,nat,sym,middle,pos,cart,mq,mq_all)

implicit none

! Local
integer :: k, counter, stat, ios, l
real :: com(cart), mmm

! Arguments
integer :: e_nz, nq, nat, cart
real :: e_zmax, e_zmin, e_dz, middle
real, allocatable :: d_charge(:), e_zmesh(:) 
real,allocatable :: qqq(:), qqq_all(:), mq(:), mq_all(:)
real, allocatable :: pos(:,:)
character*4, allocatable :: atq(:), sym(:)

! Allocate the d_charge array...
e_nz=int((e_zmax-e_zmin)/e_dz)
if (e_nz.lt.1) then
   write(99,*) "Your mesh along z (electro section) is ill-defined. Modify e_zmin, e_zmax and/or e_dz"
   stop
endif
allocate(d_charge(e_nz),e_zmesh(e_nz))
d_charge(:)=0.0d0
do k=1,e_nz
   e_zmesh(k)=e_zmin+((k*e_dz)-(e_dz/2.0d0))
enddo

! Read the charges
open(unit=975, file='charges.q', status='old')
read(975,*)
counter=0 
do while(.true.)
  read(975,*,iostat=ios)
  if( ios > 0 ) then
    stop 'problem somewhere'
  else if( ios < 0 ) then 
    exit
  else
    counter=counter+1
  end if
end do
rewind(975)
allocate(atq(counter),qqq(counter),qqq_all(nat),mq(counter),mq_all(nat))
read(975,*)
do k=1,counter
   read(975,*) atq(k), qqq(k), mq(k)
   atq(k)=adjustl(trim(atq(k)))
enddo
close(975)

nq=counter

! Pair the charges and the masses with the atoms...
do k=1,nat
   do l=1,nq
      if (sym(k).eq.atq(l)) then
         qqq_all(k)=qqq(l)
         mq_all(k)=mq(l)
         exit
      endif        
   enddo
enddo

! Get the com at the beginning, in order to scale the whole thing approximately in the middle of the box at
! every frame
! Note that because of the vacuum along z, there is no need of the fancy computation of the com
! within periodic boundary conditions, the dumb way should be enough!

com(:)=0.0
mmm=0.0
do k=1,nat
   com(cart)=com(cart)+(pos(cart,k)*mq_all(k))
   mmm=mmm+mq_all(k)
enddo
com(cart)=com(cart)/mmm
middle=com(cart)

end subroutine electro_alloc

subroutine electro(e_zmin,e_zmax,e_dz,nat,e_nz,pos,d_charge,cart,sym, &
                   atq,qqq,nq,qqq_all,middle,mq,mq_all)

implicit none

!! Arguments
integer :: nat, e_nz, cart, nq ! ns, nz, cart
!integer, allocatable :: n_ws(:), list_ws(:,:)
real :: e_zmin, e_zmax, e_dz, middle
real, allocatable :: pos(:,:), d_charge(:)
real, allocatable :: qqq(:), qqq_all(:), mq(:), mq_all(:)
character*4, allocatable :: atq(:), sym(:)
!character*3 :: switch_zdens

! Local
integer :: i, j, k
real :: lb, ub, com(cart), mmm
!character*3 :: ta
!character*1 :: fl

! Remove COM
com(:)=0.0
mmm=0.0
do j=1,nat
   com(cart)=com(cart)+(pos(cart,j)*mq_all(j))
   mmm=mmm+mq_all(j)
enddo
com(cart)=com(cart)/mmm

!! Charge density profile along z...
do j=1,nat
   ! rescale the positions wrt the center of mass of each frame - than move the whole thing in the middle
   pos(cart,j)=pos(cart,j)-com(cart)+middle
   do k=1,e_nz
      lb=e_zmin+(k-1)*e_dz
      ub=lb+e_dz
      if (pos(cart,j).gt.lb.and.pos(cart,j).le.ub) then
         d_charge(k)=d_charge(k)+qqq_all(j)
      endif
   enddo
enddo

return

end subroutine electro

end module MOD_electro

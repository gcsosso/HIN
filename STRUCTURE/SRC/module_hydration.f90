module MOD_hydration

contains

subroutine hydration_alloc(nat,ns,sym,n_ws,list_ws,o_ns,list_nw,n_nw,n_ow,o_dist,nh_bins,nh_rmax,nh_r,nh_mol,nh_atm,nh_color)

implicit none

! Arguments
integer :: nat, ns, o_ns, n_nw, n_ow, nh_bins
integer, allocatable :: n_ws(:), list_ws(:,:), list_nw(:), nh_mol(:), nh_atm(:,:), nh_color(:)
real :: nh_rmax
real, allocatable :: o_dist(:), nh_r(:)
character*4, allocatable :: sym(:)

! Local
integer :: counter, i
real :: nh_dr

do i=1,ns
  if (sym(list_ws(i,1)).eq.'OW') then
    o_ns=i
  endif
enddo

n_ow=n_ws(o_ns)
n_nw=nat-n_ws(o_ns)*4
if (.not. allocated(list_nw)) then ! If not allocated (i.e. switch_gr=no) then will need to allocate and set values
  allocate(list_nw(n_nw))
  allocate(o_dist(n_ow))
  counter=1
  do i=1,nat
    if ((sym(i).eq.'OW').or.(sym(i).eq.'HW1').or.(sym(i).eq.'HW2').or.(sym(i).eq.'MW')) then
      continue
    else
      list_nw(counter)=i
      counter=counter+1
    endif
  enddo
endif

allocate(nh_r(nh_bins),nh_mol(nh_bins),nh_atm(n_nw,nh_bins),nh_color(nat))
nh_mol(:)=0
nh_atm(:,:)=0
nh_dr=nh_rmax/(real(nh_bins))
do i=1,nh_bins
  nh_r(i)=(real(i)*nh_dr)
enddo

end subroutine hydration_alloc


subroutine hydration(resname,resnum,nat,pos,list_ws,o_ns,cart,icell,list_nw,n_nw,n_ow,o_dist,nh_bins,nh_r,nh_mol,nh_atm,nh_color)

implicit none

! Arguments
integer :: nat, o_ns, cart, n_nw, n_ow, nh_bins
integer, allocatable :: list_ws(:,:), list_nw(:), nh_mol(:), nh_atm(:,:), resnum(:), nh_color(:)
real :: icell(cart*cart)
real, allocatable :: pos(:,:), nh_r(:), o_dist(:)
character*5,allocatable :: resname(:)

! Local
integer :: i, j, i_spc, j_spc
real :: i_pos(3), j_pos(3), xdf, ydf, zdf, r_ij

do i=1,n_nw
  i_spc=list_nw(i)
  i_pos(1)=pos(1,i_spc) ; i_pos(2)=pos(2,i_spc) ; i_pos(3)=pos(3,i_spc)
  do j=1,n_ow
    j_spc=list_ws(o_ns,j)
    j_pos(1)=pos(1,j_spc) ; j_pos(2)=pos(2,j_spc) ; j_pos(3)=pos(3,j_spc)
    xdf=i_pos(1)-j_pos(1) ; ydf=i_pos(2)-j_pos(2) ; zdf=i_pos(3)-j_pos(3)

    call images(cart,0,1,1,icell,xdf,ydf,zdf)
    r_ij=sqrt(xdf**2.0d0+ydf**2.0d0+zdf**2.0d0)
    if (i.eq.1) then ! First iteration of outer loop, assign regardless of value
      o_dist(j)=r_ij
    elseif (r_ij.lt.o_dist(j)) then ! If the water is closer to another atom, assign this new shorter distance
      o_dist(j)=r_ij
    endif
  enddo
enddo

nh_color(:)=0
do i=1,n_ow
  do j=1,nh_bins
    if (o_dist(i).le.nh_r(j)) then
      nh_mol(j)=nh_mol(j)+1
      i_spc=list_ws(o_ns,i)
      nh_color(i_spc)=1
    endif
  enddo
enddo

open(unit=165, file='hin_structure.out.hydration.color', status='unknown')
write(165,*) nh_color(:)
!close(165)

end subroutine hydration


end module MOD_hydration

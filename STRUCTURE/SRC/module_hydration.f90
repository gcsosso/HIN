module MOD_hydration

contains

subroutine hydration_alloc(nat,ns,sym,n_ws,list_ws,o_ns,list_nw,n_nw,n_ow,o_dist,nh_bins,nh_rmax,nh_r,nh_mol,nh_atm,nh_color,o_nhbrs,ooo_ang,order_t)

implicit none

! Arguments
integer :: nat, ns, o_ns, n_nw, n_ow, nh_bins
integer, allocatable :: n_ws(:), list_ws(:,:), list_nw(:), nh_mol(:), nh_atm(:,:), nh_color(:), o_nhbrs(:,:)
real :: nh_rmax, ooo_ang(6)
real, allocatable :: o_dist(:), nh_r(:), order_t(:)
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

allocate(o_nhbrs(n_ow,4), order_t(n_ow))

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

! Hydration number
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

end subroutine hydration

subroutine t_order(n_ow,list_ws,o_ns,pos,cart,icell,o_nhbrs,ooo_ang,order_t,resname,resnum)

! Arguments
integer :: n_ow, o_ns, cart
integer, allocatable :: list_ws(:,:), o_nhbrs(:,:), resnum(:)
real :: ooo_ang(6), icell(cart*cart)
real, allocatable :: pos(:,:), order_t(:)
character*5,allocatable :: resname(:)

! Local
integer :: i, j, k, i_spc, j_spc, k_spc, max_loc, counter, frame
integer, parameter :: n_nn=4
real :: i_pos(3), j_pos(3), k_pos(3), xdf, ydf, zdf, r_ij, r_ik, oo_dist(4), max_dist, v_ij(3), v_ik(3), v_prod, r_prod, theta, t_sum, t_ord

open(unit=166, file='hin_structure.out.t_order', status='unknown')

! Find 4 nearest neighbours
do i=1,n_ow
  counter=1
  i_spc=list_ws(o_ns,i)
  i_pos(1)=pos(1,i_spc) ; i_pos(2)=pos(2,i_spc) ; i_pos(3)=pos(3,i_spc)

  do j=1,n_ow
    j_spc=list_ws(o_ns,j)
    if (j_spc.eq.i_spc) then
      cycle
    else
      j_pos(1)=pos(1,j_spc) ; j_pos(2)=pos(2,j_spc) ; j_pos(3)=pos(3,j_spc)
      xdf=i_pos(1)-j_pos(1) ; ydf=i_pos(2)-j_pos(2) ; zdf=i_pos(3)-j_pos(3)

      call images(cart,0,1,1,icell,xdf,ydf,zdf)
      r_ij=sqrt(xdf**2.0+ydf**2.0+zdf**2.0)
      if (counter.le.n_nn) then
        oo_dist(counter)=r_ij
        o_nhbrs(i,counter)=j_spc
        counter=counter+1
      else
        max_loc=maxloc(oo_dist,1)
        max_dist=oo_dist(max_loc)
        if (r_ij.lt.max_dist) then
          oo_dist(max_loc)=r_ij
          o_nhbrs(i,max_loc)=j_spc
        endif
      endif
    endif
  enddo
  !write(*,*) i_spc, o_nhbrs(i,:), oo_dist(:)
enddo

! Find Oi-Oj-Ok angle and calculate t parameter
do i=1,n_ow
  !counter=1
  t_sum=0.0d0
  t_ord=0.0d0
  i_spc=list_ws(o_ns,i)
  i_pos(1)=pos(1,i_spc) ; i_pos(2)=pos(2,i_spc) ; i_pos(3)=pos(3,i_spc)

  do j=1,n_nn ! Number of neighbours
    j_spc=o_nhbrs(i,j)
    j_pos(1)=pos(1,j_spc) ; j_pos(2)=pos(2,j_spc) ; j_pos(3)=pos(3,j_spc)
    xdf=j_pos(1)-i_pos(1) ; ydf=j_pos(2)-i_pos(2) ; zdf=j_pos(3)-i_pos(3)
    call images(cart,0,1,1,icell,xdf,ydf,zdf)
    v_ij(1)=xdf ; v_ij(2)=ydf ; v_ij(3)=zdf
    r_ij=sqrt(xdf**2.0+ydf**2.0+zdf**2.0)

    do k=j+1,n_nn ! Number of neighbours (k>j)
      k_spc=o_nhbrs(i,k)
      k_pos(1)=pos(1,k_spc) ; k_pos(2)=pos(2,k_spc) ; k_pos(3)=pos(3,k_spc)
      xdf=k_pos(1)-i_pos(1) ; ydf=k_pos(2)-i_pos(2) ; zdf=k_pos(3)-i_pos(3)
      call images(cart,0,1,1,icell,xdf,ydf,zdf)
      v_ik(1)=xdf ; v_ik(2)=ydf ; v_ik(3)=zdf
      r_ik=sqrt(xdf**2.0+ydf**2.0+zdf**2.0)

      v_prod=((v_ij(1)*v_ik(1))+(v_ij(2)*v_ik(2))+(v_ij(3)*v_ik(3)))
      r_prod=r_ij*r_ik
      theta=acos(v_prod/r_prod)
      t_sum=t_sum+(cos(theta)+real(1/3))**2 ! cos(acos(x))==x so could condense these lines
      !write(*,*) theta, cos(theta)
      !ooo_ang(counter)=theta
      !counter=counter+1
      !write(*,*) i_spc, j_spc, k_spc, v_prod, r_prod, counter
    enddo
  enddo
  t_ord=1.0d0-((3.0d0/8.0d0)*t_sum)
  write(166,"(i6,a6,i6,i6,i6,i6,f8.3)") resnum(i_spc), resname(i_spc), resnum(o_nhbrs(i,1)), resnum(o_nhbrs(i,2)), resnum(o_nhbrs(i,3)), resnum(o_nhbrs(i,4)), t_ord
  write(166,*)
enddo

end subroutine t_order


end module MOD_hydration

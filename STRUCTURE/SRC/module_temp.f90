module MOD_temp

contains

subroutine temp_alloc(nat) !,nh_bins,nh_rcut,nh_r,nh_mol,nh_atm,nh_color,o_nhbrs,ooo_ang,order_t,n_all_ws,list_all_ws,list_cs,n_cs)

implicit none

! Arguments
integer :: nat!, nh_bins, n_all_ws, n_cs
!integer, allocatable :: nh_mol(:), nh_atm(:,:), nh_color(:), o_nhbrs(:,:), list_all_ws(:), list_cs(:)
!real :: nh_rcut, ooo_ang(6)
!real, allocatable :: nh_r(:), order_t(:)
!
!! Local
!integer :: i
!real :: nh_dr
!
!allocate(nh_r(nh_bins),nh_mol(nh_bins),nh_atm(n_cs,nh_bins),nh_color(nat),o_nhbrs(n_all_ws,4),order_t(n_all_ws))
!nh_mol(:)=0
!nh_atm(:,:)=0
!nh_dr=nh_rcut/(real(nh_bins))
!do i=1,nh_bins
!  nh_r(i)=(real(i)*nh_dr)
!end do
!
end subroutine temp_alloc
!
!
subroutine temp(nat,pos_past,pos,lag,ts,cart,sym,zmin,zmax,dz,list_ws,n_ws)!(nh_bins,nh_r,nh_mol,nh_atm,nh_color,n_all_ws,n_filtered,list_all_ws,filt_param)

implicit none
! Arguments
real, allocatable :: pos_past(:,:), pos(:,:)
integer :: nat, lag, cart,nz
real :: ts, zmin,zmax,dz
character*4, allocatable :: sym(:)
!integer, allocatable :: nh_mol(:), nh_atm(:,:), nh_color(:), list_all_ws(:)
!real, allocatable :: nh_r(:), filt_param(:)

! Local
integer :: i, counter,l,j
real, parameter :: m_ow=2.6566962E-26 ! Kg
real, parameter :: kb=1.380649E-23 ! J K-1 
real :: vel(cart,nat), v(nat), color(nat), k(nat), t(nat), ave_T,lb,ub
integer, allocatable :: list_ws(:,:),n_ws(:)
!vel=pos-pos_past
!ave_T=0.0
!counter=0

!do i=1,nat
!  if (trim(adjustl(sym(i))).eq."OW") then
!    counter=counter+1
!    v(i)=sqrt(vel(1,i)**2+vel(2,i)**2+vel(3,i)**2)/ts ! magnitude of the velocity vector for each atom, nm / ns - or, indeed, m/s
!    k(i)=0.5*m_ow*(v(i)**2)
!    t(i)=(2.0/(3.0*kb))*k(i)
!    color(i)=t(i)
!    ave_T=ave_T+color(i)
!   else
!    color(i)=0.0
!  endif
!enddo
!
!write(*,*) "Average temperature [K] - oxygen atoms: ", ave_T/real(counter)
!open(unit=76, file='temp_color.dat', status="unknown")
!write(76,*) color(:)
!pos_past=pos

!_________z_density_____________
!!copied from Gabriele's code above don't mess with above- incorporate below into Gabriele's
!code once following checks out!
nz=int((zmax-zmin)/dz)
vel=pos-pos_past

do l=1,nz
  ave_T=0.0
  counter=0
  !

  lb=zmin+(l-1)*dz
  ub=lb+dz
   do i=1,nat
     if (trim(adjustl(sym(i))).eq."OW") then
        if (pos(3,i).gt.lb.and.pos(3,i).le.ub) then
          counter=counter+1
          v(i)=sqrt(vel(1,i)**2+vel(2,i)**2+vel(3,i)**2)/(ts) ! magnitude of the velocity vector for each atom, nm / ns - or, indeed, m/s
          k(i)=0.5*m_ow*(v(i)**2)
          !t(i)=(2.0/(3.0*kb))*k(i)
          t(i)=(3.0/(kb))*k(i) !DoF includes Settle and removal of CoM constraints
          color(i)=t(i)
          ave_T=ave_T+color(i)
        else
          color(i)=0.0
        endif
     endif
   enddo

   !
   write(*,*) "slice, avg_temp(slice)" ,l, ave_T/real(counter), ave_T,counter
enddo

open(unit=76, file='temp_color.dat', status="unknown")
write(76,*) color(:)
pos_past=pos

end subroutine temp



!
!
!subroutine t_order(pos,cart,icell,o_nhbrs,ooo_ang,order_t,t_rcut,resname,resnum,filt_max,list_filtered,n_filtered,filt_param)
!
!! Arguments
!integer :: cart, n_filtered(2)
!integer, allocatable :: o_nhbrs(:,:), resnum(:), list_filtered(:,:)
!integer, parameter :: nn=4
!real :: ooo_ang(6), icell(cart*cart), t_rcut, filt_max
!real, allocatable :: pos(:,:), order_t(:), filt_param(:)
!character*5,allocatable :: resname(:)
!
!! Local
!integer :: i, j, k, i_spc, j_spc, k_spc, max_loc, frame
!real :: i_pos(3), j_pos(3), k_pos(3), xdf, ydf, zdf, r_ij, r_ik, oo_dist(4), max_dist, v_ij(3), v_ik(3), v_prod, r_prod, theta, t_sum, t_ord
!
!open(unit=166, file='hin_structure.out.t_order', status='unknown')
!order_t(:)=0.0d0
!
!! Find 4 nearest neighbours
!do i=1,n_filtered(1)
!  oo_dist(:)=t_rcut
!  o_nhbrs(i,:)=0
!  i_spc=list_filtered(1,i)
!  i_pos(1)=pos(1,i_spc) ; i_pos(2)=pos(2,i_spc) ; i_pos(3)=pos(3,i_spc)
!
!  do j=1,n_filtered(2)
!    j_spc=list_filtered(2,j)
!
!    if (j_spc.ne.i_spc) then
!      j_pos(1)=pos(1,j_spc) ; j_pos(2)=pos(2,j_spc) ; j_pos(3)=pos(3,j_spc)
!      xdf=i_pos(1)-j_pos(1) ; ydf=i_pos(2)-j_pos(2) ; zdf=i_pos(3)-j_pos(3)
!      call images(cart,0,1,1,icell,xdf,ydf,zdf)
!      r_ij=sqrt(xdf**2.0+ydf**2.0+zdf**2.0)
!
!      max_loc=maxloc(oo_dist,1)
!      max_dist=oo_dist(max_loc)
!
!      if ((r_ij.lt.max_dist)) then
!        oo_dist(max_loc)=r_ij
!        o_nhbrs(i,max_loc)=j_spc
!      endif
!    endif
!  enddo
!  !write(*,*) i_spc, o_nhbrs(i,:), oo_dist(:)
!enddo
!
!! Find Oi-Oj-Ok angle and calculate t parameter
!do i=1,n_filtered(1)
!  !counter=1
!  t_sum=0.0d0
!  t_ord=0.0d0
!  i_spc=list_filtered(1,i)
!  i_pos(1)=pos(1,i_spc) ; i_pos(2)=pos(2,i_spc) ; i_pos(3)=pos(3,i_spc)
!
!  do j=1,nn ! Number of neighbours
!    j_spc=o_nhbrs(i,j)
!    if (j_spc.eq.0) then
!      continue
!    else
!      j_pos(1)=pos(1,j_spc) ; j_pos(2)=pos(2,j_spc) ; j_pos(3)=pos(3,j_spc)
!      xdf=j_pos(1)-i_pos(1) ; ydf=j_pos(2)-i_pos(2) ; zdf=j_pos(3)-i_pos(3)
!      call images(cart,0,1,1,icell,xdf,ydf,zdf)
!      v_ij(1)=xdf ; v_ij(2)=ydf ; v_ij(3)=zdf
!      r_ij=sqrt(xdf**2.0+ydf**2.0+zdf**2.0)
!    end if
!    do k=j+1,nn ! Number of neighbours (k>j)
!      k_spc=o_nhbrs(i,k)
!      if (j_spc.eq.0 .or. k_spc.eq.0) then
!        t_sum=t_sum+(2.0d0/3.0d0)**2
!      else
!        k_pos(1)=pos(1,k_spc) ; k_pos(2)=pos(2,k_spc) ; k_pos(3)=pos(3,k_spc)
!        xdf=k_pos(1)-i_pos(1) ; ydf=k_pos(2)-i_pos(2) ; zdf=k_pos(3)-i_pos(3)
!        call images(cart,0,1,1,icell,xdf,ydf,zdf)
!        v_ik(1)=xdf ; v_ik(2)=ydf ; v_ik(3)=zdf
!        r_ik=sqrt(xdf**2.0+ydf**2.0+zdf**2.0)
!        v_prod=((v_ij(1)*v_ik(1))+(v_ij(2)*v_ik(2))+(v_ij(3)*v_ik(3)))
!        r_prod=r_ij*r_ik
!        theta=acos(v_prod/r_prod)
!        t_sum=t_sum+(cos(theta)+(1.0d0/3.0d0))**2 ! cos(acos(x))==x so could condense these lines
!      end if
!    enddo
!  enddo
!  t_ord=1.0d0-((3.0d0/8.0d0)*t_sum)
!  order_t(i)=t_ord
!  !write(166,"(i6,a6,i6,i6,i6,i6,f8.3,f8.3)") resnum(i_spc), resname(i_spc), resnum(o_nhbrs(i,1)), resnum(o_nhbrs(i,2)), resnum(o_nhbrs(i,3)), resnum(o_nhbrs(i,4)), t_ord, filt_param(i)
!enddo
!
!write(166,*) (list_filtered(1,i), filt_param(i), order_t(i), i=1,n_filtered(1)) ! Output format: [index, distance, order]/atom for all atoms, single line per frame
!
!end subroutine t_order


end module MOD_temp

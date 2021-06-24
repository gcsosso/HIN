module MOD_hydration

contains

subroutine hydration_alloc(nat,nh_bins,nh_rcut,nh_r,nh_mol,nh_atm,nh_color,o_nhbrs,ooo_ang,order_t,n_all_ws,list_all_ws,list_cs,n_cs)

implicit none

! Arguments
integer :: nat, nh_bins, n_all_ws, n_cs
integer, allocatable :: nh_mol(:), nh_atm(:,:), nh_color(:), o_nhbrs(:,:), list_all_ws(:), list_cs(:)
real :: nh_rcut, ooo_ang(6)
real, allocatable :: nh_r(:), order_t(:)

! Local
integer :: i
real :: nh_dr

allocate(nh_r(nh_bins),nh_mol(nh_bins),nh_atm(n_cs,nh_bins),nh_color(nat),o_nhbrs(n_all_ws,4),order_t(n_all_ws))
nh_mol(:)=0
nh_atm(:,:)=0
nh_dr=nh_rcut/(real(nh_bins))
do i=1,nh_bins
  nh_r(i)=(real(i)*nh_dr)
end do

end subroutine hydration_alloc


subroutine hbond_alloc(nat,sym,resname,pos,cart,icell,ws,n_ws,hb_ws,hb_ws_filt,n_hb_x,list_hb_x,n_hb_hyd,list_hb_hyd_ws1,list_hb_hyd_ws2,sum_hb_bonds,sum_hb_filt,hb_ang)

implicit none

! Arguments
integer :: nat, cart, n_hb_x(2), sum_hb_bonds(2), sum_hb_filt
integer, allocatable :: n_ws(:), list_hb_ws(:,:), list_hb_x(:,:), n_hb_hyd(:,:), list_hb_hyd_ws1(:,:), list_hb_hyd_ws2(:,:)
real :: icell(cart*cart), hb_ang
real, allocatable :: pos(:,:)
character(4), allocatable :: sym(:), ws(:)
character(*), allocatable :: resname(:)
character(20) :: hb_ws(2)
logical(1) :: hb_ws_filt(2)

! Local
integer :: i, n_ow, n_hb_ws(2)
real, parameter :: pi=4.0d0*datan(1.0d0)
real, parameter :: deg2rad=pi/180.0d0

allocate(list_hb_ws(2,nat),list_hb_x(2,nat),n_hb_hyd(2,nat))
n_hb_ws(:)=0
hb_ws_filt(:)=.false.

if (hb_ws(1).eq.ws(1)) then ! hydration -ws1 input
  if (ws(1).ne."OW") then ! species -ws input
    write(99,*) "If using the filter with module_hydration, species -ws must be OW."
  else
    hb_ws_filt(1)=.true.
  endif
else
  call get_indices(nat,sym,resname,1,hb_ws,n_hb_ws,list_hb_ws) ! find all OW indices
  call get_hbda(nat,sym,pos,cart,icell,1,n_hb_ws,list_hb_ws,n_hb_x,list_hb_x,n_hb_hyd,list_hb_hyd_ws1) ! then find all potential Hbond donors and acceptors
endif

if (hb_ws(2).eq.hb_ws(1)) then
  if (hb_ws_filt(1)) then
    hb_ws_filt(2)=.true.
  else
    n_hb_ws(2)=n_hb_ws(1) ; n_hb_x(2)=n_hb_x(1) ; list_hb_x(2,:)=list_hb_x(1,:) ; n_hb_hyd(2,:)=n_hb_hyd(1,:) ; list_hb_hyd_ws2=list_hb_hyd_ws1
  endif
elseif (hb_ws(2).eq.ws(1)) then ! hydration -ws2 input
  if (ws(1).ne."OW") then ! species -ws input
    write(99,*) "If using the filter with module_hydration, species -ws must be OW."
  else
    hb_ws_filt(2)=.true.
  endif
else
    call get_indices(nat,sym,resname,2,hb_ws,n_hb_ws,list_hb_ws)
    call get_hbda(nat,sym,pos,cart,icell,2,n_hb_ws,list_hb_ws,n_hb_x,list_hb_x,n_hb_hyd,list_hb_hyd_ws2)
endif

! allocate(n_hb_bonds_ws1(n_hb_x(1),2))
! allocate(n_hb_bonds_ws2(n_hb_x(2),2))

! if (hb_ws_filt(1)) then
!   allocate(n_hb_bonds_ws1(n_ws(1),2))
! else
!   allocate(n_hb_bonds_ws1(n_hb_x(1),2))
! endif
! if (hb_ws_filt(2)) then
!   allocate(n_hb_bonds_ws2(n_ws(1),2))
! else
!   allocate(n_hb_bonds_ws2(n_hb_x(2),2))
! endif

deallocate(list_hb_ws)

sum_hb_bonds(:)=0
sum_hb_filt=0
hb_ang=deg2rad*hb_ang ! Convert angle from degrees to radians

open(unit=167, file='hin_structure.out.hbonds.ws1', status='unknown')
open(unit=168, file='hin_structure.out.hbonds.ws2', status='unknown')
open(unit=169, file='hin_structure.out.hbonds', status='unknown')

end subroutine hbond_alloc


subroutine get_indices(nat,sym,resname,ws_i,hb_ws,n_hb_ws,list_hb_ws)  ! Find numbers and indices from user input: atom type (e.g. OW, O1), residue name (e.g. PVA, THR) or index range (e.g. 12:12, 0:122)

  implicit none

  integer :: nat, i, ws_i, delim_index, ws_index(2), n_hb_ws(2)
  integer, allocatable :: list_hb_ws(:,:)
  character(20) :: hb_ws(2)
  character(1) :: colon=':'
  character(4), allocatable :: sym(:)
  character(*), allocatable :: resname(:)
  logical(1) :: ws_range

  if (verify(colon,hb_ws(ws_i)).eq.0) then ! If colon found then interpret as a index range
    ws_range = .true.
    delim_index = scan(hb_ws(ws_i),colon)
    read(hb_ws(ws_i)(1:delim_index-1),*) ws_index(1)
    read(hb_ws(ws_i)(delim_index+1:),*) ws_index(2)

  else ; ws_range = .false. ; end if ! Otherwise interpret as resname/atom name

  do i=1,nat
    if (trim(adjustl(hb_ws(ws_i))).eq.sym(i)) then
      n_hb_ws(ws_i)=n_hb_ws(ws_i)+1
      list_hb_ws(ws_i,n_hb_ws(ws_i))=i
    elseif (trim(adjustl(hb_ws(ws_i))).eq.resname(i)) then
      n_hb_ws(ws_i)=n_hb_ws(ws_i)+1
      list_hb_ws(ws_i,n_hb_ws(ws_i))=i
    elseif ((ws_range).and.(i.ge.ws_index(1)).and.(i.le.ws_index(2))) then ;
      n_hb_ws(ws_i)=n_hb_ws(ws_i)+1
      list_hb_ws(ws_i,n_hb_ws(ws_i))=i
    endif
  enddo

end subroutine get_indices


subroutine get_hbda(nat,sym,pos,cart,icell,ws_i,n_hb_ws,list_hb_ws,n_hb_x,list_hb_x,n_hb_hyd,list_hb_hyd) ! Filters donors/acceptors from ws inputs

  implicit none

  integer :: nat, cart, i, j, i_spc, j_spc, ws_i, n_hb_ws(2), n_hb_x(2)
  integer, allocatable :: list_hb_ws(:,:), list_hb_x(:,:), n_hb_hyd(:,:), list_hb_hyd(:,:)
  real :: icell(cart*cart), i_pos(3), j_pos(3), xdf, ydf, zdf, ij_dist
  real, allocatable :: pos(:,:)
  real, parameter :: xh_cut=0.12d0 ! Length of covalent bond between X (donor) and H
  character(4), allocatable :: sym(:)

  allocate(list_hb_hyd(n_hb_ws(ws_i),3))
  list_hb_hyd(:,:)=0
  n_hb_hyd(ws_i,:)=0
  n_hb_x(ws_i)=0

  do i=1,n_hb_ws(ws_i)
    i_spc=list_hb_ws(ws_i,i)
    if (sum(scan(sym(i_spc),["N","O","F"])).gt.0) then ! If sym(i_spc) contains N, O or F - need to test this works for different inputs
      n_hb_x(ws_i)=n_hb_x(ws_i)+1
      list_hb_x(ws_i,n_hb_x(ws_i))=i_spc
      i_pos(1)=pos(1,i_spc) ; i_pos(2)=pos(2,i_spc) ; i_pos(3)=pos(3,i_spc)
      do j_spc=1,nat
        if (scan(sym(j_spc),"H").gt.0) then
          j_pos(1)=pos(1,j_spc) ; j_pos(2)=pos(2,j_spc) ; j_pos(3)=pos(3,j_spc)
          xdf=i_pos(1)-j_pos(1) ; ydf=i_pos(2)-j_pos(2) ; zdf=i_pos(3)-j_pos(3)
          call images(cart,0,1,1,icell,xdf,ydf,zdf)
          ij_dist=sqrt(xdf**2.0+ydf**2.0+zdf**2.0)

          if (ij_dist.lt.xh_cut) then ! If distance between X and H < length of X-H covalent bond
            n_hb_hyd(ws_i,n_hb_x(ws_i))=n_hb_hyd(ws_i,n_hb_x(ws_i))+1 ! n_hb_hyd(ws_i,i) gives number of hydrogens bonded to species list_hb_x(ws_i,i)
            list_hb_hyd(n_hb_x(ws_i),n_hb_hyd(ws_i,n_hb_x(ws_i)))=j_spc ! list_hb_hyd(ws_i,j) gives index of j-th hydrogen atom bonded to list_hb_x(ws_i,i)
          endif
        endif
      enddo
    endif
  enddo

end subroutine get_hbda


subroutine hbond(sym,pos,cart,icell,hb_ws,hb_ws_filt,n_hb_x,list_hb_x,n_hb_hyd,list_hb_hyd_ws1,list_hb_hyd_ws2,n_filtered,list_filtered,sum_hb_bonds,hb_dist,hb_ang,dostuff) ! Main routine for hbonds

implicit none

! Arguments
integer :: cart, n_hb_x(2), n_filtered(2), sum_hb_bonds(2), dostuff
integer, allocatable :: list_hb_x(:,:), n_hb_hyd(:,:), list_hb_hyd_ws1(:,:), list_hb_hyd_ws2(:,:), list_filtered(:,:), n_hb_bonds_ws1(:,:), n_hb_bonds_ws2(:,:)
real :: icell(cart*cart), hb_dist, hb_ang
real, allocatable :: pos(:,:)
character(20) :: hb_ws(2)
character(4), allocatable :: sym(:)
logical(1) :: hb_ws_filt(2)

! Local
integer :: i, j, k, l, i_spc, j_spc, k_spc, l_spc, ih_spc, jh_spc, n_hb_bonds
real :: ij(3), ik(3), jk(3), il(3), jl(3), i_pos(3), j_pos(3), k_pos(3), l_pos(3), ij_dist, ik_dist, jk_dist, il_dist, jl_dist, dot_prod, cos_theta, theta

! For alanine analysis (remove later)
integer :: n_angs
real  :: xh_dist(3), xh_ang(2)
real, parameter :: pi=4.0d0*datan(1.0d0)
real, parameter :: rad2deg=180.0/pi
real, allocatable :: hb_angs(:)
allocate(hb_angs(n_filtered(1)*2))
n_angs=0
xh_dist(:)=0.0d0
hb_angs(:)=0.0d0

if (hb_ws_filt(1)) then
  allocate(n_hb_bonds_ws1(n_filtered(1),2))
else
  allocate(n_hb_bonds_ws1(n_hb_x(1),2))
endif
if (hb_ws_filt(2)) then
  allocate(n_hb_bonds_ws2(n_filtered(1),2))
else
  allocate(n_hb_bonds_ws2(n_hb_x(2),2))
endif

n_hb_bonds_ws1(:,:)=0
n_hb_bonds_ws2(:,:)=0
n_hb_bonds=0

if (hb_ws_filt(1)) then
  do i=1, n_filtered(1)
    i_spc=list_filtered(1,i)
    i_pos(1)=pos(1,i_spc) ; i_pos(2)=pos(2,i_spc) ; i_pos(3)=pos(3,i_spc)

    if (hb_ws_filt(2)) then
      do j=i+1,n_filtered(1)
        j_spc=list_filtered(1,j)
        j_pos(1)=pos(1,j_spc) ; j_pos(2)=pos(2,j_spc) ; j_pos(3)=pos(3,j_spc)
        ij(1)=i_pos(1)-j_pos(1) ; ij(2)=i_pos(2)-j_pos(2) ; ij(3)=i_pos(3)-j_pos(3)
        call images(cart,0,1,1,icell,ij(1),ij(2),ij(3))
        ij_dist=sqrt(ij(1)**2.0+ij(2)**2.0+ij(3)**2.0) ! distance between donor X and acceptor O

        if (ij_dist.lt.hb_dist) then
          do k=1,2
            ih_spc=i_spc+k ! H atom bound to Oi
            jh_spc=j_spc+k ! H atom bound to Oj
            call check_hbond(sym,pos,cart,icell,i_spc,j_spc,ih_spc,hb_ws,hb_ang,n_hb_bonds,n_hb_bonds_ws1,n_hb_bonds_ws2,i,j,1,2,theta,ik_dist,jk_dist) ! check if Oi-Hi-Oj bond meets criteria (i=donor, j=acc)
            call check_hbond(sym,pos,cart,icell,i_spc,j_spc,jh_spc,hb_ws,hb_ang,n_hb_bonds,n_hb_bonds_ws1,n_hb_bonds_ws2,i,j,2,1,theta,ik_dist,jk_dist) ! check if Oi-Hj-Oj bond meets criteria (j=donor, i=acc)
          enddo
        endif
      enddo

    else
      do j=1, n_hb_x(2)
        j_spc=list_hb_x(2,j)
        j_pos(1)=pos(1,j_spc) ; j_pos(2)=pos(2,j_spc) ; j_pos(3)=pos(3,j_spc)
        ij(1)=i_pos(1)-j_pos(1) ; ij(2)=i_pos(2)-j_pos(2) ; ij(3)=i_pos(3)-j_pos(3)
        call images(cart,0,1,1,icell,ij(1),ij(2),ij(3))
        ij_dist=sqrt(ij(1)**2.0+ij(2)**2.0+ij(3)**2.0) ! distance between donor X and acceptor O

        if (ij_dist.lt.hb_dist) then
          do k=1,2
            ih_spc=i_spc+k ! H atom bound to Oi
            call check_hbond(sym,pos,cart,icell,i_spc,j_spc,ih_spc,hb_ws,hb_ang,n_hb_bonds,n_hb_bonds_ws1,n_hb_bonds_ws2,i,j,1,2,theta,ik_dist,jk_dist) ! check if Oi-Hi-Xj bond meets criteria
          enddo
          do l=1,n_hb_hyd(2,j)
            jh_spc=list_hb_hyd_ws2(j,l) ! H atom bound to Xj
            call check_hbond(sym,pos,cart,icell,i_spc,j_spc,jh_spc,hb_ws,hb_ang,n_hb_bonds,n_hb_bonds_ws1,n_hb_bonds_ws2,i,j,2,1,theta,ik_dist,jk_dist) ! check if Oi-Hj-Xj bond meets criteria
          enddo
        endif
      enddo
    endif
  enddo

else
  do i=1, n_hb_x(1) ! Loop through potential hbond donors/acceptors
    i_spc=list_hb_x(1,i) ! i_spc=list_hb_x(2,j)
    i_pos(1)=pos(1,i_spc) ; i_pos(2)=pos(2,i_spc) ; i_pos(3)=pos(3,i_spc)

    if (hb_ws_filt(2)) then
      do j=1, n_filtered(1)
        j_spc=list_filtered(1,j)
        j_pos(1)=pos(1,j_spc) ; j_pos(2)=pos(2,j_spc) ; j_pos(3)=pos(3,j_spc)
        ij(1)=i_pos(1)-j_pos(1) ; ij(2)=i_pos(2)-j_pos(2) ; ij(3)=i_pos(3)-j_pos(3)
        call images(cart,0,1,1,icell,ij(1),ij(2),ij(3))
        ij_dist=sqrt(ij(1)**2.0+ij(2)**2.0+ij(3)**2.0) ! distance between donor X and acceptor O

        if (ij_dist.lt.hb_dist) then

          if (n_hb_hyd(1,i).gt.0) then
            n_angs=n_angs+1
            do k=1,n_hb_hyd(1,i)
              ih_spc=list_hb_hyd_ws1(i,k)! H atom bound to Xi
              call check_hbond(sym,pos,cart,icell,i_spc,j_spc,ih_spc,hb_ws,hb_ang,n_hb_bonds,n_hb_bonds_ws1,n_hb_bonds_ws2,i,j,1,2,theta,ik_dist,jk_dist) ! check if Xi-Hi-Oj bond meets criteria
              xh_dist(k)=jk_dist
              xh_ang(k)=theta
            enddo
            if ((xh_dist(1).lt.xh_dist(2)).and.(xh_dist(1).lt.xh_dist(3))) then ! Only want the angles for the closest H atom
              hb_angs(n_angs)=xh_ang(1)
            elseif ((xh_dist(2).lt.xh_dist(1)).and.(xh_dist(2).lt.xh_dist(3))) then
              hb_angs(n_angs)=xh_ang(2)
            else
              hb_angs(n_angs)=xh_ang(3)
            endif
          endif

          n_angs=n_angs+1
          do l=1,2
            jh_spc=j_spc+l ! H atom bound to Oj
            call check_hbond(sym,pos,cart,icell,i_spc,j_spc,jh_spc,hb_ws,hb_ang,n_hb_bonds,n_hb_bonds_ws1,n_hb_bonds_ws2,i,j,2,1,theta,ik_dist,jk_dist) ! check if Xi-Hj-Oj bond meets criteria
            xh_dist(l)=ik_dist
            xh_ang(l)=theta
          enddo

          if (xh_dist(1).lt.xh_dist(2)) then ! Only want the angles for the closest H atom
            hb_angs(n_angs)=xh_ang(1)
          else
            hb_angs(n_angs)=xh_ang(2)
          endif

        endif
      enddo

    else
      if (hb_ws(1).eq.(hb_ws(2))) then ! If -ws1 = -ws2 then must avoid double counting bonds
        do j=i+1, n_hb_x(2)
          j_spc=list_hb_x(2,j)
          j_pos(1)=pos(1,j_spc) ; j_pos(2)=pos(2,j_spc) ; j_pos(3)=pos(3,j_spc)
          ij(1)=i_pos(1)-j_pos(1) ; ij(2)=i_pos(2)-j_pos(2) ; ij(3)=i_pos(3)-j_pos(3)
          call images(cart,0,1,1,icell,ij(1),ij(2),ij(3))
          ij_dist=sqrt(ij(1)**2.0+ij(2)**2.0+ij(3)**2.0) ! distance between donor X and acceptor X

          if (ij_dist.lt.hb_dist) then
            do k=1,n_hb_hyd(1,i)
              ih_spc=list_hb_hyd_ws1(i,k)! H atom bound to Xi
              call check_hbond(sym,pos,cart,icell,i_spc,j_spc,ih_spc,hb_ws,hb_ang,n_hb_bonds,n_hb_bonds_ws1,n_hb_bonds_ws2,i,j,1,2,theta,ik_dist,jk_dist) ! check if Xi-Hi-Xj bond meets criteria
            enddo
            do l=1,n_hb_hyd(2,j)
              jh_spc=list_hb_hyd_ws2(j,l) ! H atom bound to Xj
              call check_hbond(sym,pos,cart,icell,i_spc,j_spc,jh_spc,hb_ws,hb_ang,n_hb_bonds,n_hb_bonds_ws1,n_hb_bonds_ws2,i,j,2,1,theta,ik_dist,jk_dist) ! check if Xi-Hj-Xj bond meets criteria
            enddo
          endif
        enddo

      else
        do j=1, n_hb_x(2)
          j_spc=list_hb_x(2,j)
          j_pos(1)=pos(1,j_spc) ; j_pos(2)=pos(2,j_spc) ; j_pos(3)=pos(3,j_spc)
          ij(1)=i_pos(1)-j_pos(1) ; ij(2)=i_pos(2)-j_pos(2) ; ij(3)=i_pos(3)-j_pos(3)
          call images(cart,0,1,1,icell,ij(1),ij(2),ij(3))
          ij_dist=sqrt(ij(1)**2.0+ij(2)**2.0+ij(3)**2.0) ! distance between donor X and acceptor X

          if (ij_dist.lt.hb_dist) then
            do k=1,n_hb_hyd(1,i)
              ih_spc=list_hb_hyd_ws1(i,k)! H atom bound to Xi
              call check_hbond(sym,pos,cart,icell,i_spc,j_spc,ih_spc,hb_ws,hb_ang,n_hb_bonds,n_hb_bonds_ws1,n_hb_bonds_ws2,i,j,1,2,theta,ik_dist,jk_dist) ! check if Xi-Hi-Xj bond meets criteria (i=donor, j=acceptor)
            enddo
            do l=1,n_hb_hyd(2,j)
              jh_spc=list_hb_hyd_ws2(j,l) ! H atom bound to Xj
              call check_hbond(sym,pos,cart,icell,i_spc,j_spc,jh_spc,hb_ws,hb_ang,n_hb_bonds,n_hb_bonds_ws1,n_hb_bonds_ws2,i,j,2,1,theta,ik_dist,jk_dist) ! check if Xi-Hj-Xj bond meets criteria (j=donor, i=acceptor)
            enddo
          endif
        enddo
      endif
    endif
  enddo
endif

! Write numbers of hbonds for each atom for every frame to hin_structure.out.hbonds.{ws1/ws2}
! Output format: [index, n_donors, n_acceptors]/atom for all potential donors or acceptors, single line per frame
if (hb_ws_filt(1)) then
  write(167,*) (list_filtered(1,i), n_hb_bonds_ws1(i,1), n_hb_bonds_ws1(i,2), i=1,n_filtered(1))
else
  write(167,*) (list_hb_x(1,i), n_hb_bonds_ws1(i,1), n_hb_bonds_ws1(i,2), i=1,n_hb_x(1))
endif

if (hb_ws_filt(2)) then
  write(168,*) (list_filtered(1,i), n_hb_bonds_ws2(i,1), n_hb_bonds_ws2(i,2), i=1,n_filtered(1))
else
  write(168,*) (list_hb_x(2,i), n_hb_bonds_ws2(i,1), n_hb_bonds_ws2(i,2), i=1,n_hb_x(2))
endif

! Output for alanine analysis
write(169,'(i8,i6,20f10.2)') dostuff, n_hb_bonds, (hb_angs(i)*rad2deg, i=1,n_angs) !frame number, num. hbonds, hbond angles

! Add frame count to running count
sum_hb_bonds(1)=sum_hb_bonds(1)+n_hb_bonds

deallocate(n_hb_bonds_ws1,n_hb_bonds_ws2)

end subroutine hbond


subroutine check_hbond(sym,pos,cart,icell,i_spc,j_spc,k_spc,hb_ws,hb_ang,n_hb_bonds,n_hb_bonds_ws1,n_hb_bonds_ws2,i,j,i_gp,j_gp,theta,ik_dist,jk_dist) ! Checks for hydrogen bond between i_spc, j_spc and k_spc using geometric criteria (hb_ang)

implicit none

integer :: cart, i_spc, j_spc, k_spc, n_hb_bonds, i, j, ac, dn, i_gp, j_gp
integer, allocatable :: n_hb_bonds_ws1(:,:), n_hb_bonds_ws2(:,:)
real :: icell(cart*cart), hb_ang, i_pos(3), j_pos(3), k_pos(3), ik(3), jk(3), ik_dist, jk_dist, dot_prod, cos_theta, theta
real, allocatable :: pos(:,:)
character(4), allocatable :: sym(:)
character(20) :: hb_ws(2)

i_pos(1)=pos(1,i_spc) ; i_pos(2)=pos(2,i_spc) ; i_pos(3)=pos(3,i_spc)
j_pos(1)=pos(1,j_spc) ; j_pos(2)=pos(2,j_spc) ; j_pos(3)=pos(3,j_spc)
k_pos(1)=pos(1,k_spc) ; k_pos(2)=pos(2,k_spc) ; k_pos(3)=pos(3,k_spc)

ik(1)=i_pos(1)-k_pos(1) ; ik(2)=i_pos(2)-k_pos(2) ; ik(3)=i_pos(3)-k_pos(3)
call images(cart,0,1,1,icell,ik(1),ik(2),ik(3))
ik_dist=sqrt(ik(1)**2.0+ik(2)**2.0+ik(3)**2.0) ! distance between donor O and hydrogen

jk(1)=j_pos(1)-k_pos(1) ; jk(2)=j_pos(2)-k_pos(2) ; jk(3)=j_pos(3)-k_pos(3)
call images(cart,0,1,1,icell,jk(1),jk(2),jk(3))
jk_dist=sqrt(jk(1)**2.0+jk(2)**2.0+jk(3)**2.0)

dot_prod=((ik(1)*jk(1))+(ik(2)*jk(2))+(ik(3)*jk(3)))
cos_theta=dot_prod/(ik_dist*jk_dist)
theta=acos(cos_theta)

if (theta.gt.hb_ang) then ! i is a H-bond donor
  n_hb_bonds_ws1(i,i_gp)=n_hb_bonds_ws1(i,i_gp)+1 ! count number of H-bond donors (i_gp=1) or acceptors (i_gp=2)
  n_hb_bonds_ws2(j,j_gp)=n_hb_bonds_ws2(j,j_gp)+1 ! count number of H-bond donors (j_gp=1) or acceptors (j_gp=2)
  n_hb_bonds=n_hb_bonds+1 ! frame count total
  if (hb_ws(1).eq.hb_ws(2)) then
    n_hb_bonds_ws1(j,j_gp)=n_hb_bonds_ws1(j,j_gp)+1
    n_hb_bonds_ws2(i,i_gp)=n_hb_bonds_ws2(i,i_gp)+1
  endif
endif

end subroutine check_hbond


subroutine h_number(nh_bins,nh_r,nh_mol,nh_atm,nh_color,n_all_ws,n_filtered,list_all_ws,filt_param)

implicit none

! Arguments
integer :: nh_bins, n_all_ws, n_filtered(:)
integer, allocatable :: nh_mol(:), nh_atm(:,:), nh_color(:), list_all_ws(:)
real, allocatable :: nh_r(:), filt_param(:)

! Local
integer :: i, j, i_spc

nh_color(:)=0
do i=1,n_filtered(1)
  do j=1,nh_bins
    if (filt_param(i).le.nh_r(j)) then
      nh_mol(j)=nh_mol(j)+1
      i_spc=list_all_ws(i)
      nh_color(i_spc)=1
    endif
  enddo
enddo

open(unit=165, file='hin_structure.out.hydration.color', status='unknown')
write(165,*) nh_color(:) ! Write color file

end subroutine h_number


subroutine t_order(pos,cart,icell,o_nhbrs,ooo_ang,order_t,t_rcut,resname,resnum,filt_max,list_filtered,n_filtered,filt_param)

! Arguments
integer :: cart, n_filtered(2)
integer, allocatable :: o_nhbrs(:,:), resnum(:), list_filtered(:,:)
integer, parameter :: nn=4
real :: ooo_ang(6), icell(cart*cart), t_rcut, filt_max
real, allocatable :: pos(:,:), order_t(:), filt_param(:)
character*5,allocatable :: resname(:)

! Local
integer :: i, j, k, i_spc, j_spc, k_spc, max_loc, frame
real :: i_pos(3), j_pos(3), k_pos(3), xdf, ydf, zdf, r_ij, r_ik, oo_dist(4), max_dist, v_ij(3), v_ik(3), v_prod, r_prod, theta, t_sum, t_ord

open(unit=166, file='hin_structure.out.t_order', status='unknown')
order_t(:)=0.0d0

! Find 4 nearest neighbours
do i=1,n_filtered(1)
  oo_dist(:)=t_rcut
  o_nhbrs(i,:)=0
  i_spc=list_filtered(1,i)
  i_pos(1)=pos(1,i_spc) ; i_pos(2)=pos(2,i_spc) ; i_pos(3)=pos(3,i_spc)

  do j=1,n_filtered(2)
    j_spc=list_filtered(2,j)

    if (j_spc.ne.i_spc) then
      j_pos(1)=pos(1,j_spc) ; j_pos(2)=pos(2,j_spc) ; j_pos(3)=pos(3,j_spc)
      xdf=i_pos(1)-j_pos(1) ; ydf=i_pos(2)-j_pos(2) ; zdf=i_pos(3)-j_pos(3)
      call images(cart,0,1,1,icell,xdf,ydf,zdf)
      r_ij=sqrt(xdf**2.0+ydf**2.0+zdf**2.0)

      max_loc=maxloc(oo_dist,1)
      max_dist=oo_dist(max_loc)

      if ((r_ij.lt.max_dist)) then
        oo_dist(max_loc)=r_ij
        o_nhbrs(i,max_loc)=j_spc
      endif
    endif
  enddo
  !write(*,*) i_spc, o_nhbrs(i,:), oo_dist(:)
enddo

! Find Oi-Oj-Ok angle and calculate t parameter
do i=1,n_filtered(1)
  !counter=1
  t_sum=0.0d0
  t_ord=0.0d0
  i_spc=list_filtered(1,i)
  i_pos(1)=pos(1,i_spc) ; i_pos(2)=pos(2,i_spc) ; i_pos(3)=pos(3,i_spc)

  do j=1,nn ! Number of neighbours
    j_spc=o_nhbrs(i,j)
    if (j_spc.eq.0) then
      continue
    else
      j_pos(1)=pos(1,j_spc) ; j_pos(2)=pos(2,j_spc) ; j_pos(3)=pos(3,j_spc)
      xdf=j_pos(1)-i_pos(1) ; ydf=j_pos(2)-i_pos(2) ; zdf=j_pos(3)-i_pos(3)
      call images(cart,0,1,1,icell,xdf,ydf,zdf)
      v_ij(1)=xdf ; v_ij(2)=ydf ; v_ij(3)=zdf
      r_ij=sqrt(xdf**2.0+ydf**2.0+zdf**2.0)
    end if
    do k=j+1,nn ! Number of neighbours (k>j)
      k_spc=o_nhbrs(i,k)
      if (j_spc.eq.0 .or. k_spc.eq.0) then
        t_sum=t_sum+(2.0d0/3.0d0)**2
      else
        k_pos(1)=pos(1,k_spc) ; k_pos(2)=pos(2,k_spc) ; k_pos(3)=pos(3,k_spc)
        xdf=k_pos(1)-i_pos(1) ; ydf=k_pos(2)-i_pos(2) ; zdf=k_pos(3)-i_pos(3)
        call images(cart,0,1,1,icell,xdf,ydf,zdf)
        v_ik(1)=xdf ; v_ik(2)=ydf ; v_ik(3)=zdf
        r_ik=sqrt(xdf**2.0+ydf**2.0+zdf**2.0)
        v_prod=((v_ij(1)*v_ik(1))+(v_ij(2)*v_ik(2))+(v_ij(3)*v_ik(3)))
        r_prod=r_ij*r_ik
        theta=acos(v_prod/r_prod)
        t_sum=t_sum+(cos(theta)+(1.0d0/3.0d0))**2 ! cos(acos(x))==x so could condense these lines
      end if
    enddo
  enddo
  t_ord=1.0d0-((3.0d0/8.0d0)*t_sum)
  order_t(i)=t_ord
  !write(166,"(i6,a6,i6,i6,i6,i6,f8.3,f8.3)") resnum(i_spc), resname(i_spc), resnum(o_nhbrs(i,1)), resnum(o_nhbrs(i,2)), resnum(o_nhbrs(i,3)), resnum(o_nhbrs(i,4)), t_ord, filt_param(i)
enddo

write(166,*) (list_filtered(1,i), filt_param(i), order_t(i), i=1,n_filtered(1)) ! Output format: [index, distance, order]/atom for all atoms, single line per frame

end subroutine t_order


end module MOD_hydration

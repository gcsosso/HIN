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

! For solute-solvent hbonds
subroutine hbond1_alloc(nat,sym,resname,pos,cart,icell,hb_ws,n_hb_x,list_hb_x,n_hb_hyd,list_hb_hyd,n_hb_bonds,list_hb_bonds,sum_hb_bonds,sum_hb_filt)

implicit none

! Arguments
integer :: nat, cart, n_hb_x, sum_hb_bonds(2), sum_hb_filt
integer, allocatable :: list_hb_ws(:), list_hb_x(:),  n_hb_hyd(:), list_hb_hyd(:,:), n_hb_bonds(:,:)
real :: icell(cart*cart)
real, allocatable :: pos(:,:)
character(4), allocatable :: sym(:)
character(*), allocatable :: resname(:)
character(20) :: hb_ws
character(20), allocatable :: list_hb_bonds(:,:)

! Local
integer :: n_hb_ws, i, j, i_spc, j_spc, delim_index, ws_index(2)
real :: i_pos(3), j_pos(3), xdf, ydf, zdf, ij_dist
real, parameter :: xh_cut=0.12d0 ! Length of covalent bond between X (donor) and H
character(1) :: colon=':'
logical(1) :: ws_range

allocate(list_hb_ws(nat)) ! list_hb_ws(1) = indices of potential donor/acceptor atoms; list_hb_ws(2) = indices of all H atoms in solute
n_hb_ws=0 ! same format as above but counts instead

! Get numbers and indices for given species. Can be given as atom type (e.g. OW, O1), residue name (e.g. PVA, THR) or index range (e.g. 12:12, 0:122)
if (verify(colon,hb_ws).eq.0) then ! If colon found then interpret as a index range
  ws_range = .true.
  delim_index = scan(hb_ws,colon)
  read(hb_ws(1:delim_index-1),*) ws_index(1)
  read(hb_ws(delim_index+1:),*) ws_index(2)
else ; ws_range = .false. ; end if ! Otherwise interpret as resname/atom name

do i=1,nat
  if (trim(adjustl(hb_ws)).eq.sym(i)) then
    n_hb_ws=n_hb_ws+1
    list_hb_ws(n_hb_ws)=i
  elseif (trim(adjustl(hb_ws)).eq.resname(i)) then
    n_hb_ws=n_hb_ws+1
    list_hb_ws(n_hb_ws)=i
  elseif ((ws_range).and.(i.ge.ws_index(1)).and.(i.le.ws_index(2))) then ;
    n_hb_ws=n_hb_ws+1
    list_hb_ws(n_hb_ws)=i
  endif
enddo

! Once selected species have been found, filter donors/acceptors
allocate(list_hb_x(n_hb_ws),n_hb_hyd(n_hb_ws),list_hb_hyd(n_hb_ws,3))
list_hb_hyd(:,:)=0
n_hb_hyd(:)=0
n_hb_x=0

do i=1,n_hb_ws
  i_spc=list_hb_ws(i)
  if (sum(scan(sym(i_spc),["N","O","F"])).gt.0) then ! If sym(i_spc) contains N, O or F - need to test this works for different inputs
    n_hb_x=n_hb_x+1
    list_hb_x(n_hb_x)=i_spc
    i_pos(1)=pos(1,i_spc) ; i_pos(2)=pos(2,i_spc) ; i_pos(3)=pos(3,i_spc)
    do j_spc=1,nat
      if (scan(sym(j_spc),"H").gt.0) then
        j_pos(1)=pos(1,j_spc) ; j_pos(2)=pos(2,j_spc) ; j_pos(3)=pos(3,j_spc)
        xdf=i_pos(1)-j_pos(1) ; ydf=i_pos(2)-j_pos(2) ; zdf=i_pos(3)-j_pos(3)
        call images(cart,0,1,1,icell,xdf,ydf,zdf)
        ij_dist=sqrt(xdf**2.0+ydf**2.0+zdf**2.0)

        if (ij_dist.lt.xh_cut) then ! If distance between X and H < length of X-H covalent bond
          n_hb_hyd(n_hb_x)=n_hb_hyd(n_hb_x)+1 ! n_hb_hyd(i) gives number of hydrogens bonded to species list_hb_x(i)
          list_hb_hyd(n_hb_x,n_hb_hyd(n_hb_x))=j_spc ! list_hb_hyd(i,j) gives index of j-th hydrogen atom bonded to list_hb_x(i)
        endif
      endif
    enddo
  endif
enddo

allocate(n_hb_bonds(n_hb_x,2),list_hb_bonds(n_hb_x,4))
n_hb_bonds(:,:)=0
sum_hb_bonds(:)=0
sum_hb_filt=0

! do i=1,n_hb_x
!     write(*,*) list_hb_x(i), sym(list_hb_x(i)), n_hb_hyd(i), list_hb_hyd(i,1), list_hb_hyd(i,2), list_hb_hyd(i,3)
! enddo

end subroutine hbond1_alloc


! For solute-solute hbonds using filter
subroutine hbond2_alloc(ns,n_ws,list_ws,sym,n_hb_x,list_hb_x,n_hb_bonds,list_hb_bonds,sum_hb_bonds,sum_hb_filt)

implicit none

! Arguments
integer :: ns, n_hb_x, sum_hb_bonds(2), sum_hb_filt
integer, allocatable :: n_ws(:), list_ws(:,:), list_hb_x(:), n_hb_bonds(:,:)
character(4), allocatable :: sym(:)
character(20), allocatable :: list_hb_bonds(:,:)

! Local
integer :: i
character(20) :: hb_ws

do i=1,ns
   if (sym(list_ws(i,1)).eq.'OW') n_hb_x=n_ws(i)
enddo

allocate(list_hb_x(n_hb_x),n_hb_bonds(n_hb_x,2),list_hb_bonds(n_hb_x,4)) ! Do we even need list_hb_x ?
n_hb_bonds(:,:)=0
sum_hb_bonds(:)=0
sum_hb_filt=0

end subroutine hbond2_alloc


! For solute-solvent hbonds
subroutine hbond1(sym,pos,cart,icell,n_hb_x,list_hb_x,n_hb_hyd,list_hb_hyd,n_filtered,list_filtered,n_hb_bonds,list_hb_bonds,sum_hb_bonds)

implicit none

! Arguments
integer :: cart, n_hb_x, n_filtered(2), sum_hb_bonds(2)
integer, allocatable :: list_hb_x(:), n_hb_hyd(:), list_hb_hyd(:,:), list_filtered(:,:), n_hb_bonds(:,:)
real :: icell(cart*cart)
real, allocatable :: pos(:,:)
character(4), allocatable :: sym(:)
character(20), allocatable :: list_hb_bonds(:,:)

! Local
integer :: i, j, k, l, i_spc, j_spc, k_spc, l_spc
real :: ij(3), ik(3), jk(3), il(3), jl(3), i_pos(3), j_pos(3), k_pos(3), l_pos(3), ij_dist, ik_dist, jk_dist, il_dist, jl_dist, dot_prod, cos_theta, theta
real, parameter :: hb_cut=0.3d0, hb_ang=2.618 ! = 150 degrees

n_hb_bonds(:,:)=0

do i=1,n_hb_x
  i_spc=list_hb_x(i)
  i_pos(1)=pos(1,i_spc) ; i_pos(2)=pos(2,i_spc) ; i_pos(3)=pos(3,i_spc)

  do j=1,n_filtered(1) ! loop over Ow atoms (filtered)
    j_spc=list_filtered(1,j)
    j_pos(1)=pos(1,j_spc) ; j_pos(2)=pos(2,j_spc) ; j_pos(3)=pos(3,j_spc)
    ij(1)=i_pos(1)-j_pos(1) ; ij(2)=i_pos(2)-j_pos(2) ; ij(3)=i_pos(3)-j_pos(3)
    call images(cart,0,1,1,icell,ij(1),ij(2),ij(3))
    ij_dist=sqrt(ij(1)**2.0+ij(2)**2.0+ij(3)**2.0) ! distance between donor X and acceptor O

    if (ij_dist.lt.hb_cut) then
      if (n_hb_hyd(i).gt.0) then
        do k=1,n_hb_hyd(i) ! loop over bonded hydrogens - check if i could be a H-bond donor
          k_spc=list_hb_hyd(i,k)
          k_pos(1)=pos(1,k_spc) ; k_pos(2)=pos(2,k_spc) ; k_pos(3)=pos(3,k_spc)
          ik(1)=i_pos(1)-k_pos(1) ; ik(2)=i_pos(2)-k_pos(2) ; ik(3)=i_pos(3)-k_pos(3)
          call images(cart,0,1,1,icell,ik(1),ik(2),ik(3))
          ik_dist=sqrt(ik(1)**2.0+ik(2)**2.0+ik(3)**2.0) ! distance between donor O and hydrogen
          jk(1)=j_pos(1)-k_pos(1) ; jk(2)=j_pos(2)-k_pos(2) ; jk(3)=j_pos(3)-k_pos(3)
          call images(cart,0,1,1,icell,jk(1),jk(2),jk(3))
          dot_prod=((ik(1)*jk(1))+(ik(2)*jk(2))+(ik(3)*jk(3)))
          jk_dist=sqrt(jk(1)**2.0+jk(2)**2.0+jk(3)**2.0)
          cos_theta=dot_prod/(ik_dist*jk_dist)
          theta=acos(cos_theta)
          if (theta.gt.hb_ang) then ! i is a H-bond donor
            n_hb_bonds(i,1)=n_hb_bonds(i,1)+1 ! count number of H-bond donors
            !list_hb_bonds(i,n_hb_bonds(i,1))=trim(char(j_spc))//','//trim(char(k_spc)) ! this will do but probably some better way of labelling H-bonds uniquely
          endif
        enddo
      endif

      do l=1,2 ! loop over water H atoms - check if i could be a H-bond acceptor
        l_spc=j_spc+l ! index HW1/HW2
        l_pos(1)=pos(1,l_spc) ; l_pos(2)=pos(2,l_spc) ; l_pos(3)=pos(3,l_spc)
        il(1)=i_pos(1)-l_pos(1) ; il(2)=i_pos(2)-l_pos(2) ; il(3)=i_pos(3)-l_pos(3)
        jl(1)=j_pos(1)-l_pos(1) ; jl(2)=j_pos(2)-l_pos(2) ; jl(3)=j_pos(3)-l_pos(3)
        call images(cart,0,1,1,icell,il(1),il(2),il(3))
        call images(cart,0,1,1,icell,jl(1),jl(2),jl(3))
        il_dist=sqrt(il(1)**2.0+il(2)**2.0+il(3)**2.0)
        jl_dist=sqrt(jl(1)**2.0+jl(2)**2.0+jl(3)**2.0)
        dot_prod=((il(1)*jl(1))+(il(2)*jl(2))+(il(3)*jl(3)))
        cos_theta=dot_prod/(il_dist*jl_dist)
        theta=acos(cos_theta)
        if (theta.gt.hb_ang) then
          n_hb_bonds(i,2)=n_hb_bonds(i,2)+1
          !list_hb_bonds(i,n_hb_bonds(i,:))=trim(char(j_spc))//','//trim(char(l_spc))
        endif
      enddo
    endif
  enddo
enddo

sum_hb_bonds(1)=sum_hb_bonds(1)+sum(n_hb_bonds(:,1))
sum_hb_bonds(2)=sum_hb_bonds(2)+sum(n_hb_bonds(:,2))

! do i=1,n_hb_x
!   write(*,*) list_hb_x(i), n_hb_bonds(i,1), n_hb_bonds(i,2)
! enddo
! write(*,*) sum(n_hb_bonds(:,:))

end subroutine hbond1


! For solute-solute hbonds using filter
subroutine hbond2(sym,pos,cart,icell,n_hb_x,list_hb_x,n_filtered,list_filtered,n_hb_bonds,list_hb_bonds,sum_hb_bonds,sum_hb_filt)

implicit none

! Arguments
integer :: cart, n_hb_x, n_filtered(2), sum_hb_bonds(2), sum_hb_filt
integer, allocatable :: list_hb_x(:), list_filtered(:,:), n_hb_bonds(:,:)
real :: icell(cart*cart)
real, allocatable :: pos(:,:)
character(4), allocatable :: sym(:)
character(20), allocatable :: list_hb_bonds(:,:)

! Local
integer :: i, j, k, l, i_spc, j_spc, k_spc, l_spc
real :: ij(3), ik(3), jk(3), il(3), jl(3), i_pos(3), j_pos(3), k_pos(3), l_pos(3), ij_dist, ik_dist, jk_dist, il_dist, jl_dist, dot_prod, cos_theta, theta
real, parameter :: hb_cut=0.3d0, hb_ang=2.70 ! = 150 degrees

n_hb_bonds(:,:)=0

do i=1,n_filtered(1)
  i_spc=list_filtered(1,i)
  i_pos(1)=pos(1,i_spc) ; i_pos(2)=pos(2,i_spc) ; i_pos(3)=pos(3,i_spc)

  do j=1,n_filtered(2)
    j_spc=list_filtered(2,j)

    if (j_spc.ne.i_spc) then
      j_pos(1)=pos(1,j_spc) ; j_pos(2)=pos(2,j_spc) ; j_pos(3)=pos(3,j_spc)
      ij(1)=i_pos(1)-j_pos(1) ; ij(2)=i_pos(2)-j_pos(2) ; ij(3)=i_pos(3)-j_pos(3)
      call images(cart,0,1,1,icell,ij(1),ij(2),ij(3))
      ij_dist=sqrt(ij(1)**2.0+ij(2)**2.0+ij(3)**2.0) ! distance between donor X and acceptor O

      if (ij_dist.lt.hb_cut) then
        do k=1,2 ! loop over hydrogens on Oi
          k_spc=i_spc+k
          k_pos(1)=pos(1,k_spc) ; k_pos(2)=pos(2,k_spc) ; k_pos(3)=pos(3,k_spc)
          ik(1)=i_pos(1)-k_pos(1) ; ik(2)=i_pos(2)-k_pos(2) ; ik(3)=i_pos(3)-k_pos(3)
          call images(cart,0,1,1,icell,ik(1),ik(2),ik(3))
          ik_dist=sqrt(ik(1)**2.0+ik(2)**2.0+ik(3)**2.0) ! distance between donor O and hydrogen
          jk(1)=j_pos(1)-k_pos(1) ; jk(2)=j_pos(2)-k_pos(2) ; jk(3)=j_pos(3)-k_pos(3)
          call images(cart,0,1,1,icell,jk(1),jk(2),jk(3))
          dot_prod=((ik(1)*jk(1))+(ik(2)*jk(2))+(ik(3)*jk(3)))
          jk_dist=sqrt(jk(1)**2.0+jk(2)**2.0+jk(3)**2.0)
          cos_theta=dot_prod/(ik_dist*jk_dist)
          theta=acos(cos_theta)
          if (theta.gt.hb_ang) then ! i is an H-bond donor
            n_hb_bonds(i,1)=n_hb_bonds(i,1)+1
            ! list_hb_bonds(i,n_hb_bonds(i,1))=trim(char(j_spc))//','//trim(char(k_spc)) ! this will do but probably some better way of labelling H-bonds uniquely
          endif
        enddo

        do l=1,2 ! loop over hydrogens on Oj
          l_spc=j_spc+l
          l_pos(1)=pos(1,l_spc) ; l_pos(2)=pos(2,l_spc) ; l_pos(3)=pos(3,l_spc)
          il(1)=i_pos(1)-l_pos(1) ; il(2)=i_pos(2)-l_pos(2) ; il(3)=i_pos(3)-l_pos(3)
          jl(1)=j_pos(1)-l_pos(1) ; jl(2)=j_pos(2)-l_pos(2) ; jl(3)=j_pos(3)-l_pos(3)
          call images(cart,0,1,1,icell,il(1),il(2),il(3))
          call images(cart,0,1,1,icell,jl(1),jl(2),jl(3))
          il_dist=sqrt(il(1)**2.0+il(2)**2.0+il(3)**2.0)
          jl_dist=sqrt(jl(1)**2.0+jl(2)**2.0+jl(3)**2.0)
          dot_prod=((il(1)*jl(1))+(il(2)*jl(2))+(il(3)*jl(3)))
          cos_theta=dot_prod/(il_dist*jl_dist)
          theta=acos(cos_theta)
          if (theta.gt.hb_ang) then
            n_hb_bonds(i,2)=n_hb_bonds(i,2)+1 ! i is an H-bond acceptor
            ! list_hb_bonds(i,n_hb_bonds(i,:))=trim(char(j_spc))//','//trim(char(l_spc))
          endif
        enddo
      endif
    endif
  enddo
enddo

sum_hb_bonds(1)=sum_hb_bonds(1)+sum(n_hb_bonds(:,1))
sum_hb_bonds(2)=sum_hb_bonds(2)+sum(n_hb_bonds(:,2))
sum_hb_filt=sum_hb_filt+n_filtered(1)

end subroutine hbond2


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

module MOD_radial

contains

subroutine radial_alloc(nat,sym,resname,cart,icell,rad_ws,rad_bins,list_rad_ws,n_rad_ws,dr,half_dr,rad,rad_norm,ws1_mol)

implicit none

! Arguments
integer :: nat, cart, rad_bins, delim_index, ws1_start, ws1_end
integer, allocatable :: list_rad_ws(:,:), n_rad_ws(:)
real :: icell(cart*cart), dr, half_dr
real, allocatable :: rad(:), rad_norm(:)
character(4), allocatable :: sym(:)
character(*), allocatable :: resname(:)
character(5) :: rad_ws(2)
logical(1) :: ws1_mol

! Local
integer :: i
real :: cell_len
character(1) :: delim=':'
logical(1) :: ws1_range

! Get numbers and indices for each group of species. First species can be given as atom type (e.g. OW, O1), residue name (e.g. PVA, THR) or index range (e.g. 12:12, 0:122). Second species must be atom type (e.g. OW)
allocate(list_rad_ws(2,nat),n_rad_ws(2))
n_rad_ws(:)=0

! -ws1 input can be provided as atom index range or residue name
if (verify(delim,rad_ws(1)).eq.0) then ! If colon detected then interpret as a index range
  ws1_range = .true.
  delim_index = scan(rad_ws(1),delim)
  read(rad_ws(1)(1:delim_index-1),*) ws1_start
  read(rad_ws(1)(delim_index+1:),*) ws1_end
else ; ws1_range = .false. ; end if ! Otherwise interpret as resname/atom name

do i=1,nat
  if (trim(adjustl(rad_ws(1))).eq.sym(i)) then
    n_rad_ws(1)=n_rad_ws(1)+1
    list_rad_ws(1,n_rad_ws(1))=i
  elseif (trim(adjustl(rad_ws(1))).eq.resname(i)) then
    n_rad_ws(1)=n_rad_ws(1)+1
    list_rad_ws(1,n_rad_ws(1))=i
    ws1_mol=.true.
  elseif ((ws1_range).and.(i.ge.ws1_start).and.(i.le.ws1_end)) then ;
    n_rad_ws(1) = n_rad_ws(1) + 1
    list_rad_ws(1,n_rad_ws(1))=i
    ws1_mol=.true.
  endif
  if (trim(adjustl(rad_ws(2))).eq.sym(i)) then
    n_rad_ws(2)=n_rad_ws(2)+1
    list_rad_ws(2,n_rad_ws(2))=i
  endif
enddo

! Build radial mesh
cell_len=icell(1)
dr=cell_len/(2.0d0*rad_bins)
half_dr=dr/2.0d0

allocate(rad(rad_bins))
do i=1,rad_bins
  rad(i)=(real(i-1)+0.5d0)*dr
enddo

allocate(rad_norm(rad_bins))
rad_norm(:)=0.0d0

end subroutine radial_alloc


subroutine radial(cart,icell,pos,rad_bins,list_rad_ws,n_rad_ws,dr,half_dr,rad,rad_norm,ws1_mol,fact)

implicit none

! Arguments
integer :: cart, rad_bins
integer, allocatable :: list_rad_ws(:,:), n_rad_ws(:)
real :: icell(cart*cart), dr, half_dr, fact
real, allocatable :: pos(:,:), rad(:), rad_norm(:)
logical(1) :: ws1_mol

! Local
integer :: i, j, k, i_spc, j_spc
integer, allocatable :: rad_sum(:)
real :: i_pos(3), j_pos(3), xdf, ydf, zdf, r_ij, r2, cell_vol, rad_tmp
real, allocatable :: rad_dist(:)
real(8), parameter :: pi=4.0d0*datan(1.0d0), pi4=4.0d0*pi

allocate(rad_sum(rad_bins),rad_dist(n_rad_ws(2)))
rad_sum(:)=0
rad_dist(:)=0.0d0
rad_tmp=0.0d0

do i=1, n_rad_ws(2) ! First loop over OW
  i_spc=list_rad_ws(2,i)
  i_pos(1)=pos(1,i_spc) ; i_pos(2)=pos(2,i_spc) ; i_pos(3)=pos(3,i_spc)

  do j=1, n_rad_ws(1) ! Second loop over atoms of interest e.g. residue
    j_spc=list_rad_ws(1,j)

    if (i_spc.ne.j_spc) then
      j_pos(1)=pos(1,j_spc) ; j_pos(2)=pos(2,j_spc) ; j_pos(3)=pos(3,j_spc)
      xdf=i_pos(1)-j_pos(1) ; ydf=i_pos(2)-j_pos(2) ; zdf=i_pos(3)-j_pos(3)
      call images(cart,0,1,1,icell,xdf,ydf,zdf)
      r_ij=sqrt(xdf**2.0d0+ydf**2.0d0+zdf**2.0d0)

      if (ws1_mol) then ! Cases where ws1 is a molecule/residue/group need to be treated differently
        if (j.eq.1) then ! First iteration of outer loop, assign regardless of value
          rad_dist(i)=r_ij
        elseif (r_ij .lt. rad_dist(i)) then ! If OW is nearer another atom, assign this new shorter distance
          rad_dist(i)=r_ij
        endif
      else ! If ws1 is a single atom then just bin
        do k=1,rad_bins
          if ((r_ij.gt.rad(k)-half_dr).and.(r_ij.le.rad(k)+half_dr)) then
            rad_sum(k)=rad_sum(k)+1
          endif
        enddo
      endif
    endif
  enddo
enddo

if (ws1_mol) then ! Binning for cases where ws1 is a molecule/residue/group of atoms
  do i=1, n_rad_ws(2)
    do j=1, rad_bins
      if ((rad_dist(i).gt.rad(j)-half_dr).and.(rad_dist(i).le.rad(j)+half_dr)) then
        rad_sum(j)=rad_sum(j)+1
      endif
    enddo
  enddo
endif

! Normalisation
cell_vol=icell(1)**3.0d0
fact=pi4*dr*(n_rad_ws(2)/cell_vol)

do i=1,rad_bins
  r2=rad(i)**2.0d0
  if (ws1_mol) then
    rad_tmp=rad_sum(i)/(fact*r2*1.0d0)
  else
    rad_tmp=rad_sum(i)/(fact*r2*n_rad_ws(1))
  endif
  rad_norm(i)=rad_norm(i)+rad_tmp ! rad_norm is passed to module_output where each bin is averaged across all frames
enddo

end subroutine radial

end module MOD_radial
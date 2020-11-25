module MOD_radial

contains

subroutine radial_alloc(nat,sym,resname,rad_ws,rad_min,rad_max,rad_bins,list_rad_ws,n_rad_ws,dr,half_dr,rad,rad_norm,ws1_mol)

implicit none

! Arguments
integer :: nat, rad_bins, delim_index, ws1_start, ws1_end
integer, allocatable :: list_rad_ws(:,:), n_rad_ws(:)
real :: dr, half_dr, rad_min, rad_max
real, allocatable :: rad(:), rad_norm(:)
character(4), allocatable :: sym(:)
character(*), allocatable :: resname(:)
character(20) :: rad_ws(2)
logical(1) :: ws1_mol

! Local
integer :: i
character(1) :: delim=':'
logical(1) :: ws1_range

! Get numbers and indices for each group of species. First species can be given as atom type (e.g. OW, O1), residue name
! (e.g. PVA, THR) or index range (e.g. 12:12, 0:122). Second species must be atom type (e.g. OW)
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
dr=(rad_max-rad_min)/rad_bins
half_dr=dr/2.0d0

allocate(rad(rad_bins))
do i=1,rad_bins
  rad(i)=((real(i-1)+0.5d0)*dr)+rad_min
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
integer :: i, j, k, i_spc, j_spc, n_images=27
integer, allocatable :: rad_sum(:)
real :: i_pos(3), j_pos(3), xdf, ydf, zdf, r_ij, r2, cell_vol, rad_tmp
real, allocatable :: rad_dist(:), min_image(:)
real(8), parameter :: pi=4.0d0*datan(1.0d0), pi4=4.0d0*pi

allocate(rad_sum(rad_bins),rad_dist(n_rad_ws(2)),min_image(n_images))
rad_sum(:)=0
rad_dist(:)=0.0d0
rad_tmp=0.0d0

! Use for finding minimum distance in minimum image/for anything without spatial extent
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

! Use for finding minimum distance in all images/for things with spatial extent
! do i=1,n_rad_ws(2)
!
!   i_spc=list_rad_ws(2,i)
!   min_image(:)=0.0d0
!   call radial_images(pos,cart,icell,i_spc,list_rad_ws,n_rad_ws(2),min_image)
!
!   do j=1,n_images
!     do k=1,rad_bins
!       if ((min_image(j).gt.rad(k)-half_dr).and.(min_image(j).le.rad(k)+half_dr)) then
!         rad_sum(k)=rad_sum(k)+1
!       endif
!     enddo
!   enddo
! enddo

! Normalisation for g(r)
! cell_vol=icell(1)**3.0d0
! fact=pi4*dr*(n_rad_ws(2)/cell_vol)
!
! do i=1,rad_bins
!   r2=rad(i)**2.0d0
!   if (ws1_mol) then
!     rad_tmp=rad_sum(i)/(fact*r2*1.0d0)
!   else
!     rad_tmp=rad_sum(i)/(fact*r2*n_rad_ws(1))
!   endif
!   rad_norm(i)=rad_norm(i)+rad_tmp ! rad_norm is passed to module_output where each bin is averaged across all frames
! enddo

! Normalisation for PDF
do i=1,rad_bins
  rad_tmp=(real(rad_sum(i))/sum(rad_sum))/dr
  rad_norm(i)=rad_norm(i)+rad_tmp
enddo

end subroutine radial

subroutine radial_images(pos,cart,icell,i_spc,j_spc_list,j_spc_count,min_image)

! Arguments
integer :: cart, i_spc, j_spc_count
integer, allocatable :: j_spc_list(:,:)
real :: icell(cart*cart)
real, allocatable :: pos(:,:), min_image(:)

! Local
integer :: x, y, z, j, image
real :: i_pos(3), j_pos(3), xdf, ydf, zdf
real, allocatable :: r_ij(:)

allocate(r_ij(j_spc_count))

i_pos(1)=pos(1,i_spc) ; i_pos(2)=pos(2,i_spc) ; i_pos(3)=pos(3,i_spc) ! position of Ow
image=1

do x=-1,1
  do y=-1,1
    do z=-1,1

      do j=1, j_spc_count ! for all atoms in protein/polymer
        j_spc=j_spc_list(1,j)
        j_pos(1)=pos(1,j_spc)+(x*icell(1)) ; j_pos(2)=pos(2,j_spc)+(y*icell(1)) ; j_pos(3)=pos(3,j_spc)+(z*icell(1))
        xdf=i_pos(1)-j_pos(1) ; ydf=i_pos(2)-j_pos(2) ; zdf=i_pos(3)-j_pos(3)
        r_ij(j)=sqrt(xdf**2.0d0+ydf**2.0d0+zdf**2.0d0)
        ! write(*,*) icell(1)*x, icell(1)*y, icell(1)*z
      enddo

      min_image(image)=minval(r_ij)
      !write(*,*) image, x, y, z, min_image(image)
      ! if (min_image(image).lt.0.2d0) then
      !   write(*,*) image, x, y, z, min_image(image)
      ! endif
      image=image+1

    enddo
  enddo
enddo

end subroutine radial_images

end module MOD_radial

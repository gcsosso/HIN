module MOD_lsi


contains

subroutine lsi_alloc(nat,sym,ns,n_ws,list_ws,o_ns,cart,icell,list_nw,n_nw,n_ow,lsi_bins,dr,half_dr,rad,lsi_mol_norm,lsi_atm_norm,o_dist)

implicit none

! Arguments

integer :: nat, ns, cart, o_ns, n_nw, n_ow, lsi_bins
integer, allocatable :: n_ws(:), list_ws(:,:), list_nw(:)
real :: icell(cart*cart), dr, half_dr
real, allocatable :: rad(:), lsi_mol_norm(:), lsi_atm_norm(:,:), o_dist(:)
character*4, allocatable :: sym(:)

! Local

integer :: i, n_oo, counter
real :: l_box


! Get list of OW indices from list_ws !! Will "OW" work for all water models??
do i=1,ns
   if (sym(list_ws(i,1)).eq.'OW') then
      o_ns=i
   endif
enddo

! Get number of non-water species to allocate list_nw. !! Hard-coded for TIP4P
! (*4 for 4 particles in model)
n_nw=nat-n_ws(o_ns)*4
allocate(list_nw(n_nw))

! Number of Ow species and O-O pairs
n_ow=n_ws(o_ns)
n_oo=((n_ow-1)**2.0+(n_ow-1))/2.0


! Get indexes of all non-water species (list_nw) !! Also hard-coded for TIP4P
! particle names
counter=1
do i=1,nat
  if ((sym(i).eq.'OW').or.(sym(i).eq.'HW1').or.(sym(i).eq.'HW2').or.(sym(i).eq.'MW')) then
    continue
  else
    list_nw(counter)=i
    counter=counter+1
  endif
enddo



! Build mesh
if (icell(1).eq.icell(5).and.icell(1).eq.icell(9)) then
l_box=icell(1)
else
write(*,*) "only works with cubic boxes"
endif
dr=l_box/(real(2.0*lsi_bins))
half_dr=dr/2.0
allocate(rad(lsi_bins))
do i=1,lsi_bins
  rad(i)=(real(i-1)+0.5)*dr
enddo
allocate(lsi_mol_norm(lsi_bins),lsi_atm_norm(n_nw,lsi_bins),o_dist(n_ow))
lsi_mol_norm(:)=0.0d0
lsi_atm_norm(:,:)=0.0d0

end subroutine lsi_alloc


subroutine lsi(pos,list_ws,o_ns,cart,icell,list_nw,n_nw,n_ow,lsi_ws,lsi_bins,dr,half_dr,rad,lsi_mol_norm,lsi_atm_norm,fact,o_dist)

implicit none


! Arguments

real :: icell(cart*cart), dr, half_dr, fact
real, allocatable :: pos(:,:), rad(:), lsi_mol_norm(:), lsi_atm_norm(:,:), o_dist(:)
integer :: o_ns, cart, n_nw, n_ow, lsi_ws, lsi_bins
integer, allocatable :: list_ws(:,:), list_nw(:)


! Local

integer :: i, j, i_spc, j_spc, ir
integer, allocatable :: lsi_mol(:,:), lsi_atm(:,:)
real :: i_pos(3), j_pos(3), xdf, ydf, zdf, r_ij
real :: r2, num_i, num_j, volume, lsi_tmp
real(8), parameter :: pi=4.0*datan(1.d0), pi4=4.0*pi

allocate(lsi_mol(n_ow,lsi_bins),lsi_atm(n_nw,lsi_bins))
lsi_mol(:,:)=0
lsi_atm(:,:)=0
lsi_tmp=0.0d0
if (lsi_ws.eq.0) then ! for O-O PCF
   do i=1,n_ow
    i_spc=list_ws(o_ns,i)
    i_pos(1)=pos(1,i_spc) ; i_pos(2)=pos(2,i_spc) ; i_pos(3)=pos(3,i_spc)
    do j=1,n_ow
      j_spc=list_ws(o_ns,j)
      if (i_spc.eq.j_spc) then
        cycle
      else
        j_pos(1)=pos(1,j_spc) ; j_pos(2)=pos(2,j_spc) ; j_pos(3)=pos(3,j_spc)
        xdf=i_pos(1)-j_pos(1)
        ydf=i_pos(2)-j_pos(2)
        zdf=i_pos(3)-j_pos(3)
        call images(cart,0,1,1,icell,xdf,ydf,zdf)
        r_ij=sqrt(xdf**2.0+ydf**2.0+zdf**2.0)
        do ir=1,lsi_bins
          if ((r_ij.gt.rad(ir)-half_dr).and.(r_ij.le.rad(ir)+half_dr)) then
            lsi_mol(i,ir)=lsi_mol(i,ir)+1
          endif
        enddo
      endif
    enddo
  enddo
endif

! Normalisation
num_j=dble(n_ow)
volume=icell(1)**3.0d0
fact=pi4*dr*(num_j/volume)
if (lsi_ws.eq.0) then ! For O-O PCF [0]
  num_i=dble(n_ow)
endif
if ((lsi_ws.eq.0)) then ! For O-O PCF [0]
  do ir=1,lsi_bins
    r2=(rad(ir))**2.0d0
    lsi_tmp=lsi_mol(i,ir)/(fact*r2*num_i)
    lsi_mol_norm(ir)=lsi_mol_norm(ir)+lsi_tmp
  enddo
endif
open(unit=134, file="i_gdr.out", status='unknown')
do i=1,lsi_bins
   write(134,*) rad(i), lsi_mol(i,ir)
enddo
stop
end subroutine lsi

end module MOD_lsi

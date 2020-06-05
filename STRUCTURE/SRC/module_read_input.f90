module MOD_read_input

contains

subroutine read_input(eflag,sfile,tfile,fframe,stride,lframe,outxtc,hw_ex,switch_zdens,ns,ws,n_ws,zmin,zmax,dz, &
                      switch_rings,rings_exe,r_zmin,r_zmax,r_ns,r_ws,n_r_ws,rcut,switch_cls,plumed_exe, &
                      switch_bonds,b_zmin,b_zmax,b_dz,b_rcut,npairs,b_bins,b_bmin,b_bmax,npairs_cn,maxr, &
                      switch_hex,switch_r_cls,r_cls_W,a_thr,maxr_RINGS,switch_cages,wcol,ohstride, &
                      vmd_exe,pmpi,cls_stat,switch_xyfes,xymin,xymax,nxy,switch_r_idx,switch_ffss,thrS, &
                      switch_electro,e_zmin,e_zmax,e_dz,switch_order,wmol,axis_1,axis_2, &
                      o_zmin,o_zmax,o_dz,switch_water,switch_hbck,hbdist,hbangle,thrSS,switch_cryo,c_rcut,nr,switch_hydration)

implicit none

integer :: stride, lframe, eflag, wcol, ohstride, pmpi, nxy
integer :: ns, r_ns, fframe, i, npairs, npairs_cn, b_bins, maxr, maxr_RINGS, nr
real :: zmin, zmax, r_zmin, r_zmax, dz, rcut, b_zmin, e_zmin, e_zmax, e_dz
real :: b_zmax, b_dz, b_bmin, b_bmax, a_thr, xymin, xymax, thrS, thrSS
real :: o_zmin, o_zmax, o_dz, hbdist, hbangle, c_rcut
real, allocatable :: b_rcut(:)
character*3 :: outxtc, hw_ex, switch_zdens, switch_hex, r_cls_W, switch_electro
character*3 :: switch_rings, switch_cls, switch_bonds, switch_r_cls, switch_order
character*3 :: switch_cages, cls_stat, switch_xyfes, switch_r_idx, switch_ffss
character*3 :: switch_water, switch_hbck, switch_cryo, switch_hydration
character*100 :: sfile, tfile, rings_exe, buffer, plumed_exe, vmd_exe
integer, allocatable, intent(out) :: n_ws(:), n_r_ws(:)
character*4 :: wmol, axis_1, axis_2
character*4, allocatable, intent(out) :: ws(:), r_ws(:)

! Read input file...
open(unit=100, file='hin_structure.in', status='old')
eflag=0
do i=1,3
   read(100,*)
enddo
! Trajectory section
read(100,*) buffer, sfile  ; if (trim(adjustl(buffer)).ne.'SFILE')  eflag=1
read(100,*) buffer, tfile  ; if (trim(adjustl(buffer)).ne.'TFILE')  eflag=1
read(100,*) buffer, fframe ; if (trim(adjustl(buffer)).ne.'FFRAME') eflag=1
read(100,*) buffer, stride ; if (trim(adjustl(buffer)).ne.'STRIDE') eflag=1
read(100,*) buffer, lframe ; if (trim(adjustl(buffer)).ne.'LFRAME') eflag=1
read(100,*) buffer, outxtc ; if (trim(adjustl(buffer)).ne.'OUTXTC') eflag=1
read(100,*) buffer, ns                ; if (trim(adjustl(buffer)).ne.'NS')     eflag=1
allocate(ws(ns),n_ws(ns)) ; n_ws(:)=0
read(100,*) buffer, (ws(i), i=1,ns)   ; if (trim(adjustl(buffer)).ne.'WS')     eflag=1
read(100,*) buffer, hw_ex  ; if (trim(adjustl(buffer)).ne.'HW_EXC') eflag=1
read(100,*) ; read(100,*)
! z-Density section
read(100,*) buffer, switch_zdens      ; if (trim(adjustl(buffer)).ne.'ZDENS')  eflag=1
read(100,*) buffer, zmin              ; if (trim(adjustl(buffer)).ne.'ZMIN')   eflag=1
read(100,*) buffer, zmax              ; if (trim(adjustl(buffer)).ne.'ZMAX')   eflag=1
read(100,*) buffer, dz                ; if (trim(adjustl(buffer)).ne.'DZ')     eflag=1
read(100,*) ; read(100,*)
! 2D xy FES section
read(100,*) buffer, switch_xyfes      ; if (trim(adjustl(buffer)).ne.'XYFES')   eflag=1
read(100,*) buffer, xymin             ; if (trim(adjustl(buffer)).ne.'XYMIN')   eflag=1
read(100,*) buffer, xymax             ; if (trim(adjustl(buffer)).ne.'XYMAX')   eflag=1
read(100,*) buffer, nxy               ; if (trim(adjustl(buffer)).ne.'NXY')     eflag=1
read(100,*) ; read(100,*); read(100,*)
! Rings section
read(100,*) buffer, switch_rings          ; if (trim(adjustl(buffer)).ne.'RINGS')  eflag=1
read(100,*) buffer, rings_exe             ; if (trim(adjustl(buffer)).ne.'R_EXE')  eflag=1
read(100,*) buffer, switch_r_idx          ; if (trim(adjustl(buffer)).ne.'R_CLS_ONLY')  eflag=1
read(100,*) buffer, r_zmin                ; if (trim(adjustl(buffer)).ne.'R_ZMIN') eflag=1
read(100,*) buffer, r_zmax                ; if (trim(adjustl(buffer)).ne.'R_ZMAX') eflag=1
read(100,*) buffer, r_ns                  ; if (trim(adjustl(buffer)).ne.'R_NS')   eflag=1
allocate(r_ws(r_ns),n_r_ws(r_ns)) ; n_r_ws(:)=0
read(100,*) buffer, (r_ws(i), i=1,r_ns)   ; if (trim(adjustl(buffer)).ne.'R_WS')   eflag=1
read(100,*) buffer, rcut                  ; if (trim(adjustl(buffer)).ne.'RCUT')   eflag=1
read(100,*) buffer, switch_hbck           ; if (trim(adjustl(buffer)).ne.'HBCK')   eflag=1
read(100,*) buffer, hbdist                ; if (trim(adjustl(buffer)).ne.'HB_DIST')   eflag=1
read(100,*) buffer, hbangle               ; if (trim(adjustl(buffer)).ne.'HB_ANGLE')  eflag=1
read(100,*) buffer, maxr                  ; if (trim(adjustl(buffer)).ne.'MAXR')   eflag=1
! Convert from actual max depth of rings search (MAXR, either 4->6 or 5->7,8,9 up to now...)
! To R.I.N.G.S. notation
if (maxr.eq.6) then
   maxr_RINGS=4
else if (maxr.eq.7.or.maxr.eq.8.or.maxr.eq.9) then
   maxr_RINGS=5
else
   write(99,*) "Values of 6,7,8 or 9 only are allowed for MAXR at the moment..."
   stop
endif
read(100,*) buffer, switch_hex            ; if (trim(adjustl(buffer)).ne.'HEX')    eflag=1
read(100,*) buffer, a_thr                 ; if (trim(adjustl(buffer)).ne.'ATHR')   eflag=1
read(100,*) buffer, switch_r_cls, r_cls_W ; if (trim(adjustl(buffer)).ne.'R_CLS')  eflag=1
read(100,*) buffer, switch_cages          ; if (trim(adjustl(buffer)).ne.'CAGES')  eflag=1
read(100,*) buffer, switch_ffss, thrS, thrSS ; if (trim(adjustl(buffer)).ne.'FFS_SURF')  eflag=1
read(100,*) buffer, wcol                  ; if (trim(adjustl(buffer)).ne.'WCOL')   eflag=1
read(100,*) ; read(100,*); read(100,*)
! Clusters section
read(100,*) buffer, switch_cls          ; if (trim(adjustl(buffer)).ne.'CLS')   eflag=1
read(100,*) buffer, plumed_exe          ; if (trim(adjustl(buffer)).ne.'C_EXE') eflag=1
read(100,*) buffer, ohstride            ; if (trim(adjustl(buffer)).ne.'OHS') eflag=1
read(100,*) buffer, vmd_exe             ; if (trim(adjustl(buffer)).ne.'VMD_EXE') eflag=1
read(100,*) buffer, pmpi                ; if (trim(adjustl(buffer)).ne.'PMPI') eflag=1
read(100,*) buffer, cls_stat            ; if (trim(adjustl(buffer)).ne.'CLS_STAT') eflag=1
read(100,*) ; read(100,*)
! Bonds section
read(100,*) buffer, switch_bonds            ; if (trim(adjustl(buffer)).ne.'BON')    eflag=1
read(100,*) buffer, b_zmin                  ; if (trim(adjustl(buffer)).ne.'B_ZMIN') eflag=1
read(100,*) buffer, b_zmax                  ; if (trim(adjustl(buffer)).ne.'B_ZMAX') eflag=1
read(100,*) buffer, b_dz                    ; if (trim(adjustl(buffer)).ne.'B_DZ')   eflag=1
npairs=factorial(ns+2-1)/(2*factorial(ns-1))
allocate(b_rcut(npairs))
if (adjustl(trim(switch_bonds)).eq."yes") then
   read(100,*) buffer, (b_rcut(i), i=1,npairs) ; if (trim(adjustl(buffer)).ne.'B_RCUT') eflag=1
else
   read(100,*) buffer
endif
read(100,*) buffer, b_bins                  ; if (trim(adjustl(buffer)).ne.'B_BINS') eflag=1
read(100,*) buffer, b_bmin                  ; if (trim(adjustl(buffer)).ne.'B_BMIN') eflag=1
read(100,*) buffer, b_bmax                  ; if (trim(adjustl(buffer)).ne.'B_BMAX') eflag=1
! If the HW1/HW2 exception id present, check that HW1 is specified as the last
! species
if (hw_ex.eq.'yes') then
   if (ws(ns).ne.'HW1') then
      write(99,*) "If the HW1/HW2 option is active, HW1 must be specified as the last atomic species!"
      stop
   endif
endif

! Electrostatic section
read(100,*) ; read(100,*)
read(100,*) buffer, switch_electro          ; if (trim(adjustl(buffer)).ne.'ELECTRO') eflag=1
read(100,*) buffer, e_zmin                  ; if (trim(adjustl(buffer)).ne.'E_ZMIN') eflag=1
read(100,*) buffer, e_zmax                  ; if (trim(adjustl(buffer)).ne.'E_ZMAX') eflag=1
read(100,*) buffer, e_dz                    ; if (trim(adjustl(buffer)).ne.'E_DZ')   eflag=1

! Order  section
read(100,*) ; read(100,*)
read(100,*) buffer, switch_order            ; if (trim(adjustl(buffer)).ne.'ORDER') eflag=1
read(100,*) buffer, switch_water            ; if (trim(adjustl(buffer)).ne.'WATER') eflag=1
read(100,*) buffer, wmol                    ; if (trim(adjustl(buffer)).ne.'WMOL') eflag=1
read(100,*) buffer, axis_1, axis_2          ; if (trim(adjustl(buffer)).ne.'AXIS') eflag=1
read(100,*) buffer, o_zmin                  ; if (trim(adjustl(buffer)).ne.'O_ZMIN') eflag=1
read(100,*) buffer, o_zmax                  ; if (trim(adjustl(buffer)).ne.'O_ZMAX') eflag=1
read(100,*) buffer, o_dz                    ; if (trim(adjustl(buffer)).ne.'O_DZ')   eflag=1

! Cryo  section
read(100,*) ; read(100,*)
read(100,*) buffer, switch_cryo            ; if (trim(adjustl(buffer)).ne.'CRYO') eflag=1
read(100,*) buffer, c_rcut                 ; if (trim(adjustl(buffer)).ne.'C_RCUT') eflag=1
read(100,*) buffer, nr                     ; if (trim(adjustl(buffer)).ne.'NR') eflag=1
read(100,*) buffer, switch_hydration       ; if (trim(adjustl(buffer)).ne.'HYDRATION') eflag=1

if (eflag.eq.1) then
   write(99,*) "Something is wrong with the input file..."
   stop
endif

return

end subroutine read_input

subroutine read_gro(sfile,nat,sym,list_ws,list_r_ws,r_color,kto,n_ws,hw_ex,switch_rings,r_ns,r_ws,n_r_ws, &
                    natformat,ns,resnum,resname,idx,dummyp,ws)

implicit none

integer :: r_ns, nat, ns, i, j, idx
integer, allocatable :: n_ws(:), n_r_ws(:), list_ws(:,:), list_r_ws(:,:), r_nper(:)
integer, allocatable :: kto(:), w_rings(:,:), r_color(:), resnum(:)
real :: dummyp
character*3 :: hw_ex, switch_zdens, switch_rings
character*5,allocatable :: resname(:)
character*100 :: sfile, natformat
character*4, allocatable :: sym(:), ws(:), r_ws(:)

! Read structure file...
open(unit=101, file=trim(adjustl(sfile)), status='old')
read(101,*)
read(101,*) nat
write(natformat,*) nat
allocate(sym(nat),list_ws(ns,nat),list_r_ws(r_ns,nat),r_color(nat),kto(nat),resname(nat),resnum(nat))
n_ws(:)=0

do i=1,nat
   read(101,'(i5,2a5,i5,3f8.3,3f8.4)') resnum(i), resname(i), sym(i), idx, dummyp, dummyp, dummyp
   sym(i)=trim(adjustl(sym(i)))
   ! HW1 = HW2 exception
   if (hw_ex.eq.'yes') then
      if (trim(adjustl(sym(i))).eq.'HW2') then
         n_ws(ns)=n_ws(ns)+1
         list_ws(ns,n_ws(ns))=i
      endif
   endif
   do j=1,ns
      if (trim(adjustl(sym(i))).eq.trim(adjustl(ws(j)))) then
         n_ws(j)=n_ws(j)+1
         list_ws(j,n_ws(j))=i
      endif
   enddo
   if (trim(adjustl(switch_rings)).eq.'yes') then
      do j=1,r_ns
         !!if (trim(adjustl(r_ws(j))).ne.'OW'.and.trim(adjustl(r_ws(j))).ne.'O3'.and.trim(adjustl(r_ws(j))).ne.'OR1'.and.trim(adjustl(r_ws(j))).ne.'OR2'.and.trim(adjustl(r_ws(j))).ne.'OR3'.and.trim(adjustl(r_ws(j))).ne.'OR4') then
         !!   write(99,*) "You'll have to implement yet another type of HB check!"
         !!   stop
         !!endif
         if (trim(adjustl(sym(i))).eq.trim(adjustl(r_ws(j)))) then
            n_r_ws(j)=n_r_ws(j)+1
            list_r_ws(j,n_r_ws(j))=i
         endif
      enddo
   endif
enddo

return

end subroutine read_gro

subroutine read_first_xtc(tfile,outxtc,xtcOfile,STAT,NATOMS,nat,xd_c,xd,xd_c_out,xd_out,STEP,time,box_trans,pos,prec,icell,cart)

use, intrinsic :: iso_c_binding, only: C_NULL_CHAR, C_PTR, c_f_pointer
use xtc

implicit none

integer :: NATOMS, STAT, STEP, nat, cart
real :: time, box_trans(cart,cart), prec, icell(cart*cart)
real, allocatable :: pos(:,:)
character*3 :: outxtc
character*100 :: tfile, xtcOfile, buffer
type(C_PTR) :: xd_c, xd_c_out
type(xdrfile), pointer :: xd, xd_out

! Set the tfile name for C.
tfile=trim(tfile)//C_NULL_CHAR
! If outxtc is true, set the xtcOfile name as well
if (trim(adjustl(outxtc)).eq.'yes') then
   xtcOfile='hin_structure.out.xtc'
   xtcOfile=trim(xtcOfile)//C_NULL_CHAR
endif

! Check sfile and tfile consistency...
STAT=read_xtc_natoms(tfile,NATOMS)
if (NATOMS.ne.nat) then
   write(*,*) NATOMS, nat
   write(99,*) "sfile does not match tfile with respect to the number of atoms!"
   stop
endif
allocate(pos(3,NATOMS))

! Open the xtc file for reading. Convert C pointer to Fortran pointer.
xd_c = xdrfile_open(tfile,"r")
call c_f_pointer(xd_c,xd)
! If outxtc is true, open the xtcOfile file as well
if (trim(adjustl(outxtc)).eq.'yes') then
   xd_c_out = xdrfile_open(xtcOfile,"w")
   call c_f_pointer(xd_c_out,xd_out)
endif

! Read first frame
STAT=read_xtc(xd,NATOMS,STEP,time,box_trans,pos,prec)
icell(1)=box_trans(1,1) ; icell(2)=box_trans(1,2) ; icell(3)=box_trans(1,3)
icell(4)=box_trans(2,1) ; icell(5)=box_trans(2,2) ; icell(6)=box_trans(2,3)
icell(7)=box_trans(3,1) ; icell(8)=box_trans(3,2) ; icell(9)=box_trans(3,3)

!
!c_11 c_12 c_13 c_21 c_22 c_23 c_31 c_32 c_33

return

end subroutine read_first_xtc

subroutine read_cls_idx(lframe,fframe,stride,C_size,C_idx,nat)

implicit none

! Arguments
integer :: lframe, fframe, stride, nat
integer, allocatable :: C_size(:), C_idx(:,:)

! Local
integer :: iostat, counter, i, maxc, j
character*100 :: buffer

maxc=lframe-fframe+1

allocate(C_size(maxc),C_idx(maxc,nat))

open(unit=69, file='idx.dat', status='old')
counter=0
do j=1,maxc
   counter=counter+1
   read(69,*) C_size(counter), (C_idx(counter,i),i=1,C_size(counter))
   !!! DEBUG
   !write(*,*) "READ", counter, C_size(counter)
   !!! DEBUG
enddo

if (counter.ne.((lframe-fframe)/stride)+1) then
   write(99,*) "Nope! idx.dat has to contain the same number of frames as the .xtc. Also, STRIDE has to be =1!"
   stop
endif

close(69)

end subroutine read_cls_idx

integer function factorial(n)

implicit none

integer, intent(in) :: n
integer :: i, ans

ans=1
do i=1,n
   ans=ans*i
enddo

factorial=ans

end function factorial

end module MOD_read_input

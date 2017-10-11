module MOD_read_input

contains

subroutine read_input(eflag,sfile,tfile,fframe,stride,lframe,hw_ex,ns,ws,n_ws,zmin,zmax,dz, &
                      it0,dt,kvec,ttau,switch_diff,switch_sisf,switch_chi4,linear,nsl,td)

implicit none
                   
integer :: stride, lframe, eflag, ns, fframe, i, it0, nsl
real :: zmin, zmax, dz, dt, kvec, ttau, linear
character*3 :: hw_ex, switch_diff, switch_sisf, switch_chi4, td
character*100 :: sfile, tfile, buffer
integer, allocatable :: n_ws(:)
character*3, allocatable :: ws(:)

! Read input file...
open(unit=100, file='hin_dynamics.in', status='old')
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
read(100,*) buffer, it0    ; if (trim(adjustl(buffer)).ne.'T_OR')   eflag=1
read(100,*) buffer, dt     ; if (trim(adjustl(buffer)).ne.'DT')     eflag=1
read(100,*) buffer, hw_ex  ; if (trim(adjustl(buffer)).ne.'HW_EXC') eflag=1
read(100,*) ; read(100,*)
! Real & reciprocal space stuff
read(100,*) buffer, td                ; if (trim(adjustl(buffer)).ne.'TD')     eflag=1
read(100,*) buffer, ns                ; if (trim(adjustl(buffer)).ne.'NS')     eflag=1      
allocate(ws(ns),n_ws(ns)) ; n_ws(:)=0
read(100,*) buffer, (ws(i), i=1,ns)   ; if (trim(adjustl(buffer)).ne.'WS')     eflag=1
read(100,*) buffer, zmin              ; if (trim(adjustl(buffer)).ne.'ZMIN')   eflag=1
read(100,*) buffer, zmax              ; if (trim(adjustl(buffer)).ne.'ZMAX')   eflag=1
read(100,*) buffer, nsl               ; if (trim(adjustl(buffer)).ne.'NSL')   eflag=1
read(100,*) buffer, dz                ; if (trim(adjustl(buffer)).ne.'DZ')     eflag=1
read(100,*) buffer, kvec              ; if (trim(adjustl(buffer)).ne.'QVEC')   eflag=1
read(100,*) buffer, ttau              ; if (trim(adjustl(buffer)).ne.'TAUT')   eflag=1
read(100,*) ; read(100,*)
! Switches
read(100,*) buffer, switch_diff         ; if (trim(adjustl(buffer)).ne.'SWITCH_DIFF')  eflag=1
read(100,*) buffer, linear              ; if (trim(adjustl(buffer)).ne.'LINEAR')       eflag=1
read(100,*) buffer, switch_sisf         ; if (trim(adjustl(buffer)).ne.'SWITCH_SISF')  eflag=1
read(100,*) buffer, switch_chi4         ; if (trim(adjustl(buffer)).ne.'SWITCH_CHI4')  eflag=1
! If the HW1/HW2 exception id present, check that HW1 is specified as the last
! species
if (hw_ex.eq.'yes') then
   if (ws(ns).ne.'HW1') then
      write(*,*) "If the HW1/HW2 option is active, HW1 must be specified as the last atomic species!"
      stop
   endif  
endif   

if (eflag.eq.1) then 
   write(*,*) "Something is wrong with the input file..."
   stop
endif

return

end subroutine read_input

subroutine read_gro(sfile,nat,sym,list_ws,n_ws,hw_ex,ns,resnum,resname,idx,dummyp,ws)

implicit none

integer :: nat, ns, i, j, resnum, idx
integer, allocatable :: n_ws(:), list_ws(:,:)
real :: dummyp
character*3 :: hw_ex
character*5 :: resname
character*100 :: sfile
character*3, allocatable :: sym(:), ws(:)

! Read structure file...
open(unit=101, file=trim(adjustl(sfile)), status='old')
read(101,*)
read(101,*) nat
allocate(sym(nat),list_ws(ns,nat))
n_ws(:)=0

do i=1,nat
   read(101,'(i5,2a5,i5,3f8.3,3f8.4)') resnum, resname, sym(i), idx, dummyp, dummyp, dummyp
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
enddo

return

end subroutine read_gro

subroutine read_first_xtc(tfile,STAT,NATOMS,nat,xd_c,xd,STEP,time,box_trans,pos,prec,cart,t0max,pos0)

use, intrinsic :: iso_c_binding, only: C_NULL_CHAR, C_PTR, c_f_pointer
use xtc

implicit none

integer :: NATOMS, STAT, STEP, nat, cart, t0max
real :: time, box_trans(cart,cart), prec
real, allocatable :: pos(:,:), pos0(:,:)
character*100 :: tfile, buffer
type(C_PTR) :: xd_c
type(xdrfile), pointer :: xd

! Set the tfile name for C.
tfile=trim(tfile)//C_NULL_CHAR

! Check sfile and tfile consistency...
STAT=read_xtc_natoms(tfile,NATOMS)
if (NATOMS.ne.nat) then
   write(*,*) "sfile does not match tfile with respect to the number of atoms!"
   stop
endif
allocate(pos(cart,NATOMS),pos0(cart,t0max*NATOMS))

! Open the xtc file for reading. Convert C pointer to Fortran pointer.
xd_c = xdrfile_open(tfile,"r")
call c_f_pointer(xd_c,xd)

! Read first frame 
STAT=read_xtc(xd,NATOMS,STEP,time,box_trans,pos,prec)

return

end subroutine read_first_xtc

end module MOD_read_input

module MOD_read_input

use MOD_read_input_args
contains

subroutine read_input(ARG_LEN, sfile, tfile, fframe, lframe, stride, switch_outxtc, switch_progress, ns, ws, &
                      switch_op, switch_q, switch_qd, switch_qt, switch_t4, switch_f, switch_th, switch_t_order, filter, centre, &
                      switch_filt_param, filt_min, filt_max, q_cut, qd_cut, qt_cut, f_cut, t_rcut, op_max_cut, max_shell, &
                      switch_rings, switch_r_col, switch_r_split, switch_hbck, switch_hex, switch_r_cls, switch_cages, &
                      switch_ffss, rings_exe, r_cls_W, r_split, r_cut, hbdist, hbangle, a_thr, thrS, thrSS, maxr, maxr_RINGS,  &
                      wcol, r_ns, r_wr, r_ws, r_wh, switch_bonds, b_dz, b_rcut, b_bmin, b_bmax, b_bins, npairs, &
                      switch_zdens, zmin, zmax, dz, switch_xyfes, xymin, xymax, nxy, &
                      switch_cls, switch_f_cls, switch_cls_stat, plumed_exe, vmd_exe, &
                      f3_imax, f3_cmax, f4_imax, f4_cmin, ohstride, pmpi, switch_electro, e_zmin, e_zmax, e_dz, &
                      switch_rad, switch_rad_cn, switch_rad_smooth, switch_rad_pdf, rad_ws, rad_bins, rad_min, rad_max, &
                      switch_nh, hb_ws, hb_dist, hb_ang, switch_temp, lag,ts,switch_solv,s_rcut)

   implicit none
   integer, parameter :: LINE_LEN=255, MAX_ARGS=31, CATEGORIES=15
   integer, intent(in) :: ARG_LEN

   ! Parsing the input file
   logical(1) :: in_arg=.false., eflag=.false., tmpflag
   character(LINE_LEN) :: buf
   integer :: io, i, j, read_loc(2), num_args(CATEGORIES)=0, num_categories=0, num_cl_args
   character(ARG_LEN) :: args(CATEGORIES, MAX_ARGS)

   ! TRAJECTORY
   character(*) :: sfile, tfile
   integer :: fframe, lframe, stride
   logical(1) :: switch_outxtc, switch_progress

   ! SPECIES
   integer :: ns
   character(4), allocatable :: ws(:)
   character(*) :: filter, centre
   real :: filt_min, filt_max
   logical(1) :: switch_filt_param

   ! ORDER
   logical(1) :: switch_op, switch_q(3:6), switch_qd(3:6), switch_qt(3:6), switch_t4, switch_f(3:4), switch_th, switch_t_order
   real :: q_cut, qd_cut, qt_cut, f_cut, t_rcut, op_max_cut
   integer :: max_shell

   ! RINGS
   logical(1) :: switch_rings, switch_r_col, switch_r_split, switch_hbck, switch_hex, switch_r_cls, switch_cages, switch_ffss
   character(*) :: rings_exe, r_cls_W
   real :: r_split, r_cut, hbdist, hbangle, a_thr, thrS, thrSS
   integer :: maxr, maxr_RINGS, wcol, r_ns
   character(5), allocatable :: r_wr(:), r_ws(:)
   integer, allocatable :: r_wh(:,:)

   ! BONDS
   logical(1) :: switch_bonds
   real :: b_dz, b_bmin, b_bmax
   real, allocatable :: b_rcut(:)
   integer :: b_bins, npairs

   ! ZDENS
   logical(1) :: switch_zdens
   real :: zmin, zmax, dz

   ! XYFES
   logical(1) :: switch_xyfes
   real :: xymin, xymax
   integer :: nxy

   ! CLUSTERS
   logical(1) :: switch_cls, switch_f_cls, switch_cls_stat
   character(*) :: plumed_exe, vmd_exe
   real :: f3_imax, f3_cmax, f4_imax, f4_cmin
   integer :: ohstride, pmpi

   ! ELECTROSTATICS
   logical(1) :: switch_electro
   real :: e_zmin, e_zmax, e_dz

   ! RADIAL
   logical(1) :: switch_rad, switch_rad_cn, switch_rad_smooth, switch_rad_pdf
   character(20) :: rad_ws(2)
   integer :: rad_bins
   real :: rad_min, rad_max

   ! HYDRATION
   logical(1) :: switch_nh
   real :: hb_dist, hb_ang
   character(20) :: hb_ws(2)

   ! TEMP
   logical(1) :: switch_temp
   integer :: lag
   real :: ts

   ! SOLVATION
   logical(1) :: switch_solv
   real :: s_rcut

   read_loc(:) = 0

   num_cl_args = COMMAND_ARGUMENT_COUNT()

   if (num_cl_args.ne.0) then
      do i=1,num_cl_args ; CALL GET_COMMAND_ARGUMENT(i,args(i,1)) ; end do
      if (trim(adjustl(args(1,1))).eq.'order') then
         switch_op = .true.
         do i=2,num_cl_args
            tmpflag = .false.
            call read_traj_arg(args(i,1), tmpflag, .false._1, sfile, tfile, fframe, lframe, stride, &
                               switch_outxtc, switch_progress)
            if (.not.tmpflag) cycle
            tmpflag = .false.
            call read_ws_arg(args(i,j), eflag, filter, filt_min, filt_max, ns, ws, centre)
            call read_order_arg(args(i,1), tmpflag, .false._1, switch_q, switch_qd, switch_qt, switch_t4, switch_f, &
                                switch_th, switch_t_order, q_cut, qd_cut, qt_cut, f_cut, t_rcut, max_shell)
            if (tmpflag) then ; eflag = .true.
               write(99,*) "I don't understand the command line argument: "//trim(args(i,1))
            end if
         end do
      end if
      if (eflag) then ; write(99,*) "Something is wrong with the command line arguments..." ; stop
      else ; return ; end if
   end if

   ! If command line arguments aren't provided, we try to read input file...
   open(unit=100, file='hin_structure.in', status='old')

   do
      read(100, '(A)', iostat=io) buf ; if (io.ne.0) exit
      do i=1,LINE_LEN
         if (buf(i:i).eq.'#') exit ! Ignore anything after a hash (comments)

         ! This is a line that is not being ignored
         if (buf(i:i).eq.' ') then
            in_arg = .false.
         else if (in_arg) then
            args(read_loc(1),read_loc(2)) = trim(args(read_loc(1),read_loc(2)))//buf(i:i)
         else
            if (buf(i:i).eq.'-') then
               read_loc(2) = read_loc(2) + 1
            else
               read_loc(1) = read_loc(1) + 1
               read_loc(2) = 1
            end if
            num_categories = read_loc(1)
            num_args(read_loc(1)) = read_loc(2)
            in_arg = .true.
            args(read_loc(1),read_loc(2)) = buf(i:i)
         end if
      end do
   end do
   close(100)

   allocate(ws(MAX_ARGS)) ; ns = 0

   do i=1,num_categories
      if (args(i,1).eq.'') exit
      if (args(i,1).eq.'trajectory') then
         do j=2,num_args(i) ; if (args(i,j).eq.'') exit
            call read_traj_arg(args(i,j), eflag, .true._1, sfile, tfile, fframe, lframe, stride, &
                               switch_outxtc, switch_progress)
         end do
      else if (args(i,1).eq.'species') then
         do j=2,num_args(i) ; if (args(i,j).eq.'') exit
            call read_ws_arg(args(i,j), eflag, filter, filt_min, filt_max, ns, ws, centre)
         end do
      else if (args(i,1).eq.'order') then
         switch_op = .true.
         do j=2,num_args(i) ; if (args(i,j).eq.'') exit
            call read_order_arg(args(i,j), eflag, .true._1, switch_q, switch_qd, switch_qt, switch_t4, switch_f, &
                                switch_th, switch_t_order, q_cut, qd_cut, qt_cut, f_cut, t_rcut, max_shell)
         end do
      else if (args(i,1).eq.'rings') then
         switch_rings = .true.
         do j=2,num_args(i) ; if (args(i,j).eq.'') exit
            call read_rings_arg(args(i,j), eflag, .true._1, switch_r_col, rings_exe, r_cls_W, &
                                switch_r_split, switch_hbck, switch_hex, switch_r_cls, switch_cages, switch_ffss, &
                                r_split, r_cut, hbdist, hbangle, a_thr, thrS, thrSS, maxr, wcol)
         end do
         call read_rings_input(eflag, r_ns, r_wr, r_ws, r_wh, maxr, maxr_RINGS)
      else if (args(i,1).eq.'bonds') then
         switch_bonds = .true.
         npairs=factorial(ns+2-1)/(2*factorial(ns-1))
         allocate(b_rcut(npairs)) ; b_rcut(:) = 0.0
         j = 2 ;  do while (j.le.num_args(i)) ; if (args(i,j).eq.'') exit
            call read_bonds_arg(args(i,j), eflag, .true._1, b_dz, b_rcut, b_bins, b_bmin, b_bmax)
            if (b_rcut(1).eq.-1.0) then
               call read_b_rcut(eflag, args(i,j+1:j+npairs), b_rcut, npairs)
               j = j+npairs
            end if
            j = j+1
         end do
      else if (args(i,1).eq.'zdens') then
         switch_zdens = .true.
         do j=2,num_args(i) ; if (args(i,j).eq.'') exit
            call read_zdens_arg(args(i,j), eflag, .true._1, zmin, zmax, dz)
         end do
      else if (args(i,1).eq.'xyfes') then
         switch_xyfes = .true.
         do j=2,num_args(i) ; if (args(i,j).eq.'') exit
            call read_xyfes_arg(args(i,j), eflag, .true._1, xymin, xymax, nxy)
         end do
      else if (args(i,1).eq.'clusters') then
         do j=2,num_args(i) ; if (args(i,j).eq.'') exit
            call read_clusters_arg(args(i,j), eflag, .true._1, switch_cls, switch_f_cls, switch_f, switch_cls_stat, &
                                   plumed_exe, vmd_exe, f3_imax, f3_cmax, f4_imax, f4_cmin, ohstride, pmpi)
         end do
      else if (args(i,1).eq.'electro') then
         switch_electro = .true.
         do j=2,num_args(i) ; if (args(i,j).eq.'') exit
            call read_electro_arg(args(i,j), eflag, .true._1, e_zmin, e_zmax, e_dz)
         end do
      else if (args(i,1).eq.'radial') then
         switch_rad = .true.
         do j=2,num_args(i) ; if (args(i,j).eq.'') exit
            call read_radial_arg(args(i,j), eflag, .true._1, rad_ws, rad_bins, rad_min, rad_max, switch_rad_cn, switch_rad_smooth, switch_rad_pdf)
         end do
      else if (args(i,1).eq.'hydration') then
         do j=2,num_args(i) ; if (args(i,j).eq.'') exit
            switch_nh = .true.
            call read_hydration_arg(args(i,j), eflag, .true._1, hb_ws, hb_dist, hb_ang)
         end do
      else if (args(i,1).eq.'temp') then
         switch_temp = .true.
         do j=2,num_args(i) ; if (args(i,j).eq.'') exit
            call read_temp_arg(args(i,j), eflag, .true._1, lag, ts)
         end do
      else if (args(i,1).eq.'solv') then !defines your instruction 'solv'
         switch_solv = .true.
         do j=2,num_args(i) ; if (args(i,j).eq.'') exit
            call read_solv_arg(args(i,j),eflag, .true._1 ,s_rcut)
         end do
      else ; eflag = .true. ; write(99,*) "I don't understand the argument: "//trim(args(i,1)) ; end if
   end do

   call set_op_max_cut(switch_qd, switch_qt, q_cut, qd_cut, qt_cut, f_cut, t_rcut, op_max_cut)

   if (ns.eq.0) then
      ns = 2
      ws(1) = 'OW' ; ws(2) = 'HW'
   end if

   if (filter.eq.'z'.or.filter.eq.'shell') switch_filt_param = .true.

   if (eflag) then ; write(99,*) "Something is wrong with the input file..." ; stop ; end if

end subroutine read_input


subroutine read_gro(sfile,nat,sym,list_ws,list_r_ws,r_color,kto,switch_rings,r_ns,r_ws,r_wr,n_r_ws, &
                    natformat,ns,resnum,resname,idx,dummyp,ws,list_f_ow,n_f_ow,switch_op,coloring,list_s_ws,current_coord)

   implicit none

   integer :: r_ns, nat, ns, i, j, idx, n_f_ow, tmp_ws_len
   integer, allocatable :: n_r_ws(:), list_ws(:,:),list_r_ws(:,:),list_s_ws(:,:)
   integer, allocatable :: list_f_ow(:),n_ws(:)
   integer, allocatable :: kto(:), r_color(:),resnum(:),coloring(:),current_coord(:)
   character(5) :: tmp_ws
   real :: dummyp
   logical(1) ::switch_rings, switch_op, tmp_ws_ast
   character(5), allocatable :: resname(:)
   character(100) :: sfile, natformat
   character(4), allocatable :: sym(:), ws(:)
   character(5), allocatable :: r_ws(:), r_wr(:)

   ! Read structure file...
   open(unit=101, file=trim(adjustl(sfile)), status='old')
   read(101,*)
   read(101,*) nat
   write(natformat,*) nat
   allocate(sym(nat),list_r_ws(r_ns,nat),r_color(nat),kto(nat),resname(nat),resnum(nat),n_r_ws(r_ns))
   allocate(list_f_ow(nat))
   allocate(coloring(nat))
   allocate(n_ws(ns))
   allocate(list_s_ws(ns,nat))
   allocate(current_coord(2))

   n_r_ws(:) = 0
   n_f_ow = 0
   n_ws(:) = 0

   do i=1,nat
      read(101,'(i5,2a5,i5,3f8.3,3f8.4)') resnum(i), resname(i), sym(i), idx, dummyp, dummyp, dummyp
      sym(i)=trim(adjustl(sym(i)))
      if (switch_rings) then
         do j=1,r_ns
            !!if (trim(adjustl(r_ws(j))).ne.'OW'.and.trim(adjustl(r_ws(j))).ne.'O3'.and.trim(adjustl(r_ws(j))).ne.'OR1'.and. &
               !!trim(adjustl(r_ws(j))).ne.'OR2'.and.trim(adjustl(r_ws(j))).ne.'OR3'.and.trim(adjustl(r_ws(j))).ne.'OR4') then
            !!   write(99,*) "You'll have to implement yet another type of HB check!"
            !!   stop
            !!endif
            if ((trim(adjustl(sym(i))).eq.trim(adjustl(r_ws(j)))) &
                .and.(trim(adjustl(resname(i))).eq.trim(adjustl(r_wr(j))))) then
               n_r_ws(j)=n_r_ws(j)+1
               list_r_ws(j,n_r_ws(j))=i
            endif
         enddo
      endif

      do j=1,ns
         if (trim(adjustl(sym(i))).eq.trim(adjustl(ws(j)))) then
            n_ws(j)=n_ws(j)+1
            list_s_ws(j,n_ws(j))=i
         endif
      enddo

   enddo

   close(101)

   return

end subroutine read_gro

subroutine read_first_xtc(tfile, switch_outxtc, xtcOfile, STAT, NATOMS, nat, xd_c, xd, xd_c_out, xd_out, STEP, time, &
                          box_trans, pos, prec, icell, cart)

   use, intrinsic :: iso_c_binding, only: C_NULL_CHAR, C_PTR, c_f_pointer
   use xtc

   implicit none

   integer :: NATOMS, STAT, STEP, nat, cart
   real :: time, box_trans(cart,cart), prec, icell(cart*cart)
   real, allocatable :: pos(:,:)
   logical(1) :: switch_outxtc
   character*100 :: tfile, xtcOfile
   type(C_PTR) :: xd_c, xd_c_out
   type(xdrfile), pointer :: xd, xd_out

   ! Set the tfile name for C.
   tfile=trim(tfile)//C_NULL_CHAR
   ! If switch_outxtc is true, set the xtcOfile name as well
   if (switch_outxtc) then
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
   ! If switch_outxtc is true, open the xtcOfile file as well
   if (switch_outxtc) then
      xd_c_out = xdrfile_open(xtcOfile,"w")
      call c_f_pointer(xd_c_out,xd_out)
   endif

   ! Read first frame
   STAT=read_xtc(xd,NATOMS,STEP,time,box_trans,pos,prec)
   call box_trans2icell(cart,box_trans,icell)
!   icell(1)=box_trans(1,1) ; icell(2)=box_trans(1,2) ; icell(3)=box_trans(1,3)
!   icell(4)=box_trans(2,1) ; icell(5)=box_trans(2,2) ; icell(6)=box_trans(2,3)
!   icell(7)=box_trans(3,1) ; icell(8)=box_trans(3,2) ; icell(9)=box_trans(3,3)

   return

end subroutine read_first_xtc

subroutine box_trans2icell(cart,box_trans,icell)

   implicit none

   integer :: cart
   real :: box_trans(cart,cart), icell(cart*cart)

   icell(1)=box_trans(1,1) ; icell(2)=box_trans(1,2) ; icell(3)=box_trans(1,3)
   icell(4)=box_trans(2,1) ; icell(5)=box_trans(2,2) ; icell(6)=box_trans(2,3)
   icell(7)=box_trans(3,1) ; icell(8)=box_trans(3,2) ; icell(9)=box_trans(3,3)

   return

end subroutine box_trans2icell


subroutine read_cls_idx(lframe,fframe,stride,C_size,C_idx,nat)

   implicit none

   ! Arguments
   integer :: lframe, fframe, stride, nat
   integer, allocatable :: C_size(:), C_idx(:,:)

   ! Local
   integer :: counter, i, maxc, j

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

module MOD_read_input

contains

subroutine read_input(ARG_LEN, sfile, tfile, fframe, lframe, stride, switch_outxtc, &
                      switch_op, switch_q, switch_qd, switch_qt, switch_t4, switch_f, switch_th, &
                      filter, filt_min, filt_max, q_cut, qd_cut, qt_cut, f_cut, max_shell, op_species, &
                      switch_rings, switch_r_split, switch_hbck, switch_hex, switch_r_cls, switch_cages, switch_ffss, &
                      rings_exe, r_cls_W, r_split, r_cut, hbdist, hbangle, a_thr, thrS, thrSS, maxr, maxr_RINGS, wcol, &
                      r_ns, r_wr, r_ws, r_wh)

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
   logical(1) :: switch_outxtc
   
   ! ORDER
   logical(1) :: switch_op, switch_q(3:6), switch_qd(3:6), switch_qt(3:6), switch_t4, switch_f(3:4), switch_th
   character(*) :: filter, op_species
   real :: filt_min, filt_max, q_cut, qd_cut, qt_cut, f_cut
   integer :: max_shell
   
   ! RINGS
   logical(1) :: switch_rings, switch_r_split, switch_hbck, switch_hex, switch_r_cls, switch_cages, switch_ffss
   character(*) :: rings_exe, r_cls_W
   real :: r_split, r_cut, hbdist, hbangle, a_thr, thrS, thrSS
   integer :: maxr, maxr_RINGS, wcol, r_ns
   character(5), allocatable :: r_wr(:), r_ws(:)
   integer, allocatable :: r_wh(:,:)
   
   read_loc(:) = 0
   
   num_cl_args = COMMAND_ARGUMENT_COUNT()

   if (num_cl_args.ne.0) then
      do i=1,num_cl_args ; CALL GET_COMMAND_ARGUMENT(i,args(1,i)) ; end do
      if (trim(adjustl(args(1,1))).eq.'order') then
         switch_op = .true.
         do i=2,num_cl_args
            tmpflag = .false.
            call read_traj_arg(args(1,i), tmpflag, .false._1, sfile, tfile, fframe, lframe, stride, switch_outxtc, &
                               filter,filt_min,filt_max)
            if (.not.tmpflag) cycle
            tmpflag = .false.
            call read_order_arg(args(1,i), tmpflag, .false._1, switch_q, switch_qd, switch_qt, switch_t4, switch_f, & 
                                switch_th, q_cut, qd_cut, qt_cut, f_cut, max_shell, op_species)
            if (tmpflag) then ; eflag = .true.
               write(99,*) "I don't understand the command line argument: "//trim(args(1,i))
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
   
   do i=1,num_categories
      if (args(i,1).eq.'') exit
      if (args(i,1).eq.'trajectory') then
         do j=2,num_args(i) ; if (args(i,j).eq.'') exit
            call read_traj_arg(args(i,j), eflag, .true._1, sfile, tfile, fframe, lframe, stride, switch_outxtc, &
                               filter, filt_min, filt_max)
         end do
      else if (args(i,1).eq.'order') then
         switch_op = .true.
         do j=2,num_args(i) ; if (args(i,j).eq.'') exit
            call read_order_arg(args(i,j), eflag, .true._1, switch_q, switch_qd, switch_qt, switch_t4, switch_f, &
                                switch_th, q_cut, qd_cut, qt_cut, f_cut, max_shell, op_species)
         end do
      else if (args(i,1).eq.'rings') then
         switch_rings = .true.
         do j=2,num_args(i) ; if (args(i,j).eq.'') exit
            call read_rings_arg(args(i,j), .true._1, eflag, rings_exe, r_cls_W, r_ns, r_wr, r_ws, r_wh, &
                                switch_r_split, switch_hbck, switch_hex, switch_r_cls, switch_cages, switch_ffss, &
                                r_split, r_cut, hbdist, hbangle, a_thr, thrS, thrSS, maxr, maxr_RINGS, wcol)
         end do
      else ; eflag = .true. ; write(99,*) "I don't understand the argument: "//trim(args(i,1)) ; end if
   end do
   
   if (eflag) then ; write(99,*) "Something is wrong with the input file..." ; stop ; end if

end subroutine read_input


subroutine read_gro(sfile,nat,sym,list_ws,list_r_ws,r_color,kto,n_ws,switch_hw_ex,switch_rings,r_ns,r_ws,r_wr,n_r_ws, &
                    natformat,ns,resnum,resname,idx,dummyp,ws,list_f_ow,n_f_ow,switch_op)

   implicit none

   integer :: r_ns, nat, ns, i, j, idx, n_f_ow
   integer, allocatable :: n_ws(:), n_r_ws(:), list_ws(:,:), list_r_ws(:,:)
   integer, allocatable :: list_f_ow(:)
   integer, allocatable :: kto(:), r_color(:), resnum(:)
   real :: dummyp
   logical(1) ::switch_rings, switch_op, switch_hw_ex
   character(5),allocatable :: resname(:)
   character(100) :: sfile, natformat
   character(4), allocatable :: sym(:), ws(:)
   character(5), allocatable :: r_ws(:), r_wr(:)
   
   ! Read structure file...
   open(unit=101, file=trim(adjustl(sfile)), status='old')
   read(101,*)
   read(101,*) nat
   write(natformat,*) nat
   allocate(sym(nat),list_ws(ns,nat),list_r_ws(r_ns,nat),r_color(nat),kto(nat),resname(nat),resnum(nat))
   allocate(list_f_ow(nat))
   n_ws(:)=0
   n_f_ow = 0
   
   do i=1,nat
      read(101,'(i5,2a5,i5,3f8.3,3f8.4)') resnum(i), resname(i), sym(i), idx, dummyp, dummyp, dummyp
      sym(i)=trim(adjustl(sym(i)))
      ! HW1 = HW2 exception
      if (switch_hw_ex) then
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
      if (switch_rings) then
         do j=1,r_ns
            !!if (trim(adjustl(r_ws(j))).ne.'OW'.and.trim(adjustl(r_ws(j))).ne.'O3'.and.trim(adjustl(r_ws(j))).ne.'OR1'.and.trim(adjustl(r_ws(j))).ne.'OR2'.and.trim(adjustl(r_ws(j))).ne.'OR3'.and.trim(adjustl(r_ws(j))).ne.'OR4') then
            !!   write(99,*) "You'll have to implement yet another type of HB check!"
            !!   stop
            !!endif
            if ((trim(adjustl(sym(i))).eq.trim(adjustl(r_ws(j)))).and.(trim(adjustl(resname(i))).eq.trim(adjustl(r_wr(j))))) then
               n_r_ws(j)=n_r_ws(j)+1
               list_r_ws(j,n_r_ws(j))=i
            endif
         enddo 
      endif  
      if (switch_op) then
        if (trim(adjustl(sym(i))).eq.'OW') then
           n_f_ow=n_f_ow+1
           list_f_ow(n_f_ow)=i
        endif
      endif  
   enddo
   
   close(101)

   return

end subroutine read_gro

subroutine read_first_xtc(tfile,switch_outxtc,xtcOfile,STAT,NATOMS,nat,xd_c,xd,xd_c_out,xd_out,STEP,time,box_trans,pos,prec,icell,cart)

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
   icell(1)=box_trans(1,1) ; icell(2)=box_trans(1,2) ; icell(3)=box_trans(1,3)
   icell(4)=box_trans(2,1) ; icell(5)=box_trans(2,2) ; icell(6)=box_trans(2,3)
   icell(7)=box_trans(3,1) ; icell(8)=box_trans(3,2) ; icell(9)=box_trans(3,3)
   
   return

end subroutine read_first_xtc

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

subroutine read_traj_arg(arg,eflag,log_errors,sfile,tfile,fframe,lframe,stride,switch_outxtc,filter,filt_min,filt_max)
   
   implicit none
   
   logical(1) :: eflag, switch_outxtc, log_errors
   character(*) :: arg, sfile, tfile, filter
   integer :: fframe, lframe, stride
   real :: filt_min, filt_max
   
   if (arg(1:7).eq.'-sfile=') then ; call read_arg(arg(8:), 0, 0.0, sfile, 'str', 'sfile', eflag)
   else if (arg(1:7).eq.'-tfile=') then ; call read_arg(arg(8:), 0, 0.0, tfile, 'str', 'tfile', eflag)
   else if (arg(1:8).eq.'-fframe=') then ; call read_arg(arg(9:), fframe, 0.0, '', 'int', 'fframe', eflag)
   else if (arg(1:8).eq.'-lframe=') then ; call read_arg(arg(9:), lframe, 0.0, '', 'int', 'lframe', eflag)
   else if (arg(1:8).eq.'-stride=') then ; call read_arg(arg(9:), stride, 0.0, '', 'int', 'stride', eflag)
   else if (trim(adjustl(arg)).eq.'--noxtc') then ; switch_outxtc = .false.
   else if (arg(1:8).eq.'-filter=') then ; call read_arg(arg(9:), 0, 0.0, filter, 'str', 'filter', eflag)
   else if (arg(1:5).eq.'-min=') then ; call read_arg(arg(6:), 0, filt_min, '', 'real', 'min', eflag)
   else if (arg(1:5).eq.'-max=') then ; call read_arg(arg(6:), 0, filt_max, '', 'real', 'max', eflag)
   else ; eflag = .true. ; if (log_errors) write(99,*) "I don't understand the argument: trajectory "//trim(arg) ; end if
   
end subroutine read_traj_arg

subroutine read_order_arg(arg, eflag, log_errors, switch_q, switch_qd, switch_qt, switch_t4, switch_f, &
                          switch_th, q_cut, qd_cut, qt_cut, f_cut, max_shell, op_species)
   
   implicit none
   
   logical(1) :: eflag, log_errors, switch_q(3:6), switch_qd(3:6), switch_qt(3:6), switch_t4, switch_f(3:4), switch_th
   character(*) :: arg, op_species
   real :: q_cut, qd_cut, qt_cut, f_cut
   integer :: max_shell
   
   if (trim(adjustl(arg)).eq.'--q3') then ; switch_q(3) = .true.
   else if (trim(adjustl(arg)).eq.'--q4') then ; switch_q(4) = .true.
   else if (trim(adjustl(arg)).eq.'--q6') then ; switch_q(6) = .true.
   else if (trim(adjustl(arg)).eq.'--q3d') then ; switch_qd(3) = .true.
   else if (trim(adjustl(arg)).eq.'--q4d') then ; switch_qd(4) = .true.
   else if (trim(adjustl(arg)).eq.'--q6d') then ; switch_qd(6) = .true. ; if (switch_qd(3)) switch_t4=.true.
   else if (trim(adjustl(arg)).eq.'--q3t') then ; switch_qt(3) = .true.
   else if (trim(adjustl(arg)).eq.'--q4t') then ; switch_qt(4) = .true.
   else if (trim(adjustl(arg)).eq.'--q6t') then ; switch_qt(6) = .true.
   else if (trim(adjustl(arg)).eq.'--f3') then ; switch_f(3) = .true.
   else if (trim(adjustl(arg)).eq.'--f4') then ; switch_f(4) = .true.
   else if (trim(adjustl(arg)).eq.'--th') then ; switch_th = .true.
   else if (arg(1:6).eq.'-qcut=') then ; call read_arg(arg(7:), 0, q_cut, '', 'real', 'qcut', eflag)
   else if (arg(1:7).eq.'-qdcut=') then ; call read_arg(arg(8:), 0, qd_cut, '', 'real', 'qdcut', eflag)
   else if (arg(1:7).eq.'-qtcut=') then ; call read_arg(arg(8:), 0, qt_cut, '', 'real', 'qtcut', eflag)
   else if (arg(1:6).eq.'-fcut=') then ; call read_arg(arg(7:), 0, f_cut, '', 'real', 'fcut', eflag)
   else if (arg(1:10).eq.'-maxshell=') then ; call read_arg(arg(11:), max_shell, 0.0, '', 'int', 'maxshell', eflag)
   else if (arg(1:9).eq.'-species=') then ; call read_arg(arg(10:), 0, 0.0, op_species, 'str', 'species', eflag)
   else ; eflag = .true. ; if (log_errors) write(99,*) "I don't understand the argument: order "//trim(arg) ; end if
   
end subroutine read_order_arg

subroutine read_rings_arg(arg, eflag, log_errors, rings_exe, r_cls_W, r_ns, r_wr, r_ws, r_wh, &
                          switch_r_split, switch_hbck, switch_hex, switch_r_cls, switch_cages, switch_ffss, &
                          r_split, r_cut, hbdist, hbangle, a_thr, thrS, thrSS, maxr, maxr_RINGS, wcol)
   
   implicit none
   
   logical(1) :: in_arg
   character(31) :: buf
   character(5) :: current_arg
   integer :: io, i, j, line_loc
   
   logical(1) :: eflag, log_errors, switch_r_split, switch_hbck, switch_hex, switch_r_cls, switch_cages, switch_ffss
   character(*) :: arg, rings_exe, r_cls_W
   real :: r_split, r_cut, hbdist, hbangle, a_thr, thrS, thrSS
   integer :: maxr, maxr_RINGS, wcol, r_ns
   character(5), allocatable :: r_wr(:), r_ws(:)
   integer, allocatable :: r_wh(:,:)
   
   if (arg(1:5).eq.'-exe=') then ; call read_arg(arg(6:), 0, 0.0, rings_exe, 'str', 'exe', eflag)
   else if (arg(1:7).eq.'-split=') then; switch_r_split = .true. ; call read_arg(arg(8:), 0, r_split, '', 'real', 'split', eflag)
   else if (arg(1:5).eq.'-cut=') then ; call read_arg(arg(6:), 0, r_cut, '', 'real', 'cut', eflag)
   else if (trim(adjustl(arg)).eq.'--hbck') then ; switch_hbck = .true.
   else if (arg(1:6).eq.'-dcut=') then ; call read_arg(arg(7:), 0, hbdist, '', 'real', 'dcut', eflag)
   else if (arg(1:6).eq.'-acut=') then ; call read_arg(arg(7:), 0, hbangle, '', 'real', 'acut', eflag)
   else if (arg(1:6).eq.'-maxr=') then ; call read_arg(arg(7:), maxr, 0.0, '', 'int', 'maxr', eflag)
   else if (trim(adjustl(arg)).eq.'--hex') then ; switch_hex = .true.
   else if (arg(1:7).eq.'-a_thr=') then ; call read_arg(arg(8:), 0, a_thr, '', 'real', 'a_thr', eflag)
   else if (arg(1:5).eq.'-cls=') then ; switch_r_cls = .true. ; call read_arg(arg(6:), 0, 0.0, r_cls_W, 'str', 'cls', eflag)
   else if (trim(adjustl(arg)).eq.'--cages') then ; switch_cages = .true.
   else if (trim(adjustl(arg)).eq.'--ffs_surf') then ; switch_ffss = .true.
   else if (arg(1:6).eq.'-thrS=') then ; call read_arg(arg(7:), 0, thrS, '', 'real', 'thrS', eflag)
   else if (arg(1:7).eq.'-thrSS=') then ; call read_arg(arg(8:), 0, thrSS, '', 'real', 'thrSS', eflag)
   else if (arg(1:6).eq.'-wcol=') then ; call read_arg(arg(7:), wcol, 0.0, '', 'int', 'wcol', eflag)
   else ; eflag = .true. ; if (log_errors) write(99,*) "I don't understand the argument: rings "//trim(arg) ; end if
   
   ! Convert from actual max depth of rings search (MAXR, either 4->6 or 5->7,8,9 up to now...)
   ! To R.I.N.G.S. notation
   if (maxr.eq.6) then ; maxr_RINGS=4 ; else if (maxr.eq.7.or.maxr.eq.8.or.maxr.eq.9) then ; maxr_RINGS=5
   else ; write(99,*) "Values of 6,7,8 or 9 only are allowed for MAXR at the moment..." ; stop ; endif
   
   ! Read rings input file...
   open(unit=100, file='hin_structure.rings.in', status='old')
   do ; read(100, '(A)', iostat=io) buf ; if (io.ne.0) exit ; r_ns = r_ns+1 ; end do ; close(100)
   
   allocate(r_wr(r_ns), r_ws(r_ns), r_wh(r_ns, 4)) ; r_wh = 0
   
   open(unit=100, file='hin_structure.rings.in', status='old')
   do i=1,r_ns ; read(100, '(A)', iostat=io) buf ; if (io.ne.0) exit ; line_loc = 0 ; in_arg=.false.
      do j=1,31 ; if (buf(j:j).eq.'#') then ; exit ! Ignore anything after a hash (comments)
         else 
            ! This is a line that is not being ignored
            if (buf(j:j).eq.' ') then
               if (in_arg) then
                  in_arg = .false.
                  if (line_loc.eq.1) then ; r_wr(i) = current_arg
                  else if (line_loc.eq.2) then ; r_ws(i) = current_arg
                  else if (line_loc.le.6) then
                     read(current_arg, *, iostat=io) r_wh(i, line_loc-2)
                     if (io.ne.0) then ; eflag=.true. ; exit ; end if
                  else ; exit ; end if
               end if
            else if (in_arg) then
               current_arg = trim(current_arg)//buf(j:j)
            else
               in_arg = .true.
               current_arg = buf(j:j)
               line_loc = line_loc + 1
            end if
         end if
      end do
   end do ; close(100)
   
end subroutine read_rings_arg


subroutine read_arg(arg, n, x, s, argtype, argname, eflag)
   implicit none

   logical(1) :: eflag
   character(*) :: arg, argname, argtype, s
   integer :: n, io
   real :: x

   if (argtype.eq.'int') then ; read(arg, *, iostat=io) n
   else if (argtype.eq.'real') then ; read(arg, *, iostat=io) x
   else if (argtype.eq.'str') then ; read(arg, *, iostat=io) s ; end if
   if (io.ne.0) then
      print *, trim(arg)//" is not a valid value for "//trim(argname)//"!"
      eflag=.true.
   end if

end subroutine read_arg


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


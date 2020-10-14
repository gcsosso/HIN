module MOD_read_input_args

contains

subroutine read_traj_arg(arg, eflag, log_errors, sfile, tfile, fframe, lframe, stride, &
                         switch_outxtc, switch_progress)

   implicit none

   logical(1) :: eflag, log_errors, switch_outxtc, switch_progress
   character(*) :: arg, sfile, tfile
   integer :: fframe, lframe, stride

   if (arg(1:7).eq.'-sfile=') then ; call read_arg(arg(8:), 0, 0.0, sfile, 'str', 'sfile', eflag)
   else if (arg(1:7).eq.'-tfile=') then ; call read_arg(arg(8:), 0, 0.0, tfile, 'str', 'tfile', eflag)
   else if (arg(1:8).eq.'-fframe=') then ; call read_arg(arg(9:), fframe, 0.0, '', 'int', 'fframe', eflag)
   else if (arg(1:8).eq.'-lframe=') then ; call read_arg(arg(9:), lframe, 0.0, '', 'int', 'lframe', eflag)
   else if (arg(1:8).eq.'-stride=') then ; call read_arg(arg(9:), stride, 0.0, '', 'int', 'stride', eflag)
   else if (trim(adjustl(arg)).eq.'--noxtc') then ; switch_outxtc = .false.
   else if (trim(adjustl(arg)).eq.'--progress') then ; switch_progress = .true.
   else ; eflag = .true. ; if (log_errors) write(99,*) "I don't understand the argument: trajectory "//trim(arg) ; end if

end subroutine read_traj_arg


subroutine read_ws_arg(arg, eflag, filter, filt_min, filt_max, ns, ws, centre)
   implicit none

   character(*) :: arg, filter, centre
   logical(1) :: eflag
   integer :: ns
   real :: filt_min, filt_max
   character(*), allocatable :: ws(:)

   if (arg(1:8).eq.'-filter=') then ; call read_arg(arg(9:), 0, 0.0, filter, 'str', 'filter', eflag)
   else if (arg(1:5).eq.'-min=') then ; call read_arg(arg(6:), 0, filt_min, '', 'real', 'min', eflag)
   else if (arg(1:5).eq.'-max=') then ; call read_arg(arg(6:), 0, filt_max, '', 'real', 'max', eflag)
   else if (arg(1:8).eq.'-centre=') then ; call read_arg(arg(9:), 0, 0.0, centre, 'str', 'centre', eflag)
   else ; ns = ns+1 ; call read_arg(arg(2:), 0, 0.0, ws(ns), 'str', 'species', eflag)
   end if

end subroutine read_ws_arg

subroutine read_order_arg(arg, eflag, log_errors, switch_q, switch_qd, switch_qt, switch_t4, switch_f, &
                          switch_th, switch_t_order, q_cut, qd_cut, qt_cut, f_cut, t_rcut, max_shell)

   implicit none

   logical(1) :: eflag, log_errors
   logical(1) :: switch_q(3:6), switch_qd(3:6), switch_qt(3:6), switch_t4, switch_f(3:4), switch_th, switch_t_order
   character(*) :: arg
   real :: q_cut, qd_cut, qt_cut, f_cut, t_rcut
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
   else if (trim(adjustl(arg)).eq.'--t') then ; switch_t_order = .true.
   else if (arg(1:6).eq.'-qcut=') then ; call read_arg(arg(7:), 0, q_cut, '', 'real', 'qcut', eflag)
   else if (arg(1:7).eq.'-qdcut=') then ; call read_arg(arg(8:), 0, qd_cut, '', 'real', 'qdcut', eflag)
   else if (arg(1:7).eq.'-qtcut=') then ; call read_arg(arg(8:), 0, qt_cut, '', 'real', 'qtcut', eflag)
   else if (arg(1:6).eq.'-fcut=') then ; call read_arg(arg(7:), 0, f_cut, '', 'real', 'fcut', eflag)
   else if (arg(1:6).eq.'-tcut=') then ; call read_arg(arg(7:), 0, t_rcut, '', 'real', 'tcut', eflag)
   else if (arg(1:10).eq.'-maxshell=') then ; call read_arg(arg(11:), max_shell, 0.0, '', 'int', 'maxshell', eflag)
   else ; eflag = .true. ; if (log_errors) write(99,*) "I don't understand the argument: order "//trim(arg) ; end if

end subroutine read_order_arg

subroutine read_rings_arg(arg, eflag, log_errors, switch_r_col, rings_exe, r_cls_W, &
                          switch_r_split, switch_hbck, switch_hex, switch_r_cls, switch_cages, switch_ffss, &
                          r_split, r_cut, hbdist, hbangle, a_thr, thrS, thrSS, maxr, wcol)

   implicit none

   logical(1) :: eflag, log_errors
   logical(1) :: switch_r_col, switch_r_split, switch_hbck, switch_hex, switch_r_cls, switch_cages, switch_ffss
   character(*) :: arg, rings_exe, r_cls_W
   real :: r_split, r_cut, hbdist, hbangle, a_thr, thrS, thrSS
   integer :: maxr, wcol

   if (arg(1:5).eq.'-exe=') then ; call read_arg(arg(6:), 0, 0.0, rings_exe, 'str', 'exe', eflag)
   else if (arg(1:7).eq.'-split=') then; switch_r_split = .true. ; call read_arg(arg(8:), 0, r_split, '', 'real', 'split', eflag)
   else if (arg(1:5).eq.'-cut=') then ; call read_arg(arg(6:), 0, r_cut, '', 'real', 'cut', eflag)
   else if (trim(adjustl(arg)).eq.'--hbck') then ; switch_hbck = .true.
   else if (trim(adjustl(arg)).eq.'--color') then ; switch_r_col = .true.
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

end subroutine read_rings_arg


subroutine read_bonds_arg(arg, eflag, log_errors, b_dz, b_rcut, b_bins, b_bmin, b_bmax)

   implicit none

   logical(1) :: eflag, log_errors
   character(*) :: arg
   real :: b_dz, b_bmin, b_bmax
   real, allocatable :: b_rcut(:)
   integer :: b_bins

   if (arg(1:4).eq.'-dz=') then ; call read_arg(arg(5:), 0, b_dz, '', 'real', 'dz', eflag)
   else if (trim(adjustl(arg)).eq.'-rcut') then ; b_rcut(1) = -1.0
   else if (arg(1:6).eq.'-bins=') then ; call read_arg(arg(7:), b_bins, 0.0, '', 'int', 'bins', eflag)
   else if (arg(1:6).eq.'-bmin=') then ; call read_arg(arg(7:), 0, b_bmin, '', 'real', 'bmin', eflag)
   else if (arg(1:6).eq.'-bmax=') then ; call read_arg(arg(7:), 0, b_bmax, '', 'real', 'bmax', eflag)
   else ; eflag = .true. ; if (log_errors) write(99,*) "I don't understand the argument: bonds "//trim(arg) ; end if

end subroutine read_bonds_arg


subroutine read_zdens_arg(arg, eflag, log_errors, zmin, zmax, dz)

   implicit none

   logical(1) :: eflag, log_errors
   character(*) :: arg
   real :: zmin, zmax, dz

   if (arg(1:4).eq.'-dz=') then ; call read_arg(arg(5:), 0, dz, '', 'real', 'dz', eflag)
   else if (arg(1:5).eq.'-min=') then ; call read_arg(arg(6:), 0, zmin, '', 'real', 'min', eflag)
   else if (arg(1:5).eq.'-max=') then ; call read_arg(arg(6:), 0, zmax, '', 'real', 'max', eflag)
   else ; eflag = .true. ; if (log_errors) write(99,*) "I don't understand the argument: zdens "//trim(arg) ; end if

end subroutine read_zdens_arg


subroutine read_xyfes_arg(arg, eflag, log_errors, xymin, xymax, nxy)

   implicit none

   logical(1) :: eflag, log_errors
   character(*) :: arg
   real :: xymin, xymax
   integer :: nxy

   if (arg(1:5).eq.'-nxy=') then ; call read_arg(arg(6:), nxy, 0.0, '', 'int', 'dz', eflag)
   else if (arg(1:7).eq.'-xymin=') then ; call read_arg(arg(8:), 0, xymin, '', 'real', 'xymin', eflag)
   else if (arg(1:7).eq.'-xymax=') then ; call read_arg(arg(8:), 0, xymax, '', 'real', 'xymax', eflag)
   else ; eflag = .true. ; if (log_errors) write(99,*) "I don't understand the argument: xyfes "//trim(arg) ; end if

end subroutine read_xyfes_arg


subroutine read_clusters_arg(arg, eflag, log_errors, switch_cls, switch_f_cls, switch_f, switch_cls_stat, &
                             plumed_exe, vmd_exe, f3_imax, f3_cmax, f4_imax, f4_cmin, ohstride, pmpi)

   implicit none

   logical(1) :: eflag, log_errors, switch_cls, switch_f_cls, switch_f(3:4), switch_cls_stat
   character(*) :: arg, plumed_exe, vmd_exe
   real :: f3_imax, f3_cmax, f4_imax, f4_cmin
   integer :: ohstride, pmpi

   if (trim(adjustl(arg)).eq.'--ice') then ; switch_cls = .true.
   else if (trim(adjustl(arg)).eq.'--clathrate') then ; switch_f_cls = .true. ; switch_f(3) = .true. ; switch_f(4) = .true.
   else if (trim(adjustl(arg)).eq.'--cls_stat') then ; switch_cls_stat = .true.
   else if (arg(1:12).eq.'-plumed_exe=') then ; call read_arg(arg(13:), 0, 0.0, plumed_exe, 'str', 'plumed_exe', eflag)
   else if (arg(1:9).eq.'-vmd_exe=') then ; call read_arg(arg(10:), 0, 0.0, vmd_exe, 'str', 'vmd_exe', eflag)
   else if (arg(1:11).eq.'-ice_f3max=') then ; call read_arg(arg(12:), 0, f3_imax, '', 'real', 'ice_f3max', eflag)
   else if (arg(1:17).eq.'-clathrate_f3max=') then ; call read_arg(arg(18:), 0, f3_cmax, '', 'real', 'clathrate_f3max', eflag)
   else if (arg(1:11).eq.'-ice_f4max=') then ; call read_arg(arg(12:), 0, f4_imax, '', 'real', 'ice_f4max', eflag)
   else if (arg(1:17).eq.'-clathrate_f4min=') then ; call read_arg(arg(18:), 0, f4_cmin, '', 'real', 'clathrate_f4min', eflag)
   else if (arg(1:10).eq.'-ohstride=') then ; call read_arg(arg(11:), ohstride, 0.0, '', 'int', 'ohstride', eflag)
   else if (arg(1:6).eq.'-pmpi=') then ; call read_arg(arg(7:), pmpi, 0.0, '', 'int', 'pmpi', eflag)
   else ; eflag = .true. ; if (log_errors) write(99,*) "I don't understand the argument: clusters "//trim(arg) ; end if

end subroutine read_clusters_arg


subroutine read_electro_arg(arg, eflag, log_errors, e_zmin, e_zmax, e_dz)

   implicit none

   logical(1) :: eflag, log_errors
   character(*) :: arg
   real :: e_zmin, e_zmax, e_dz

   if (arg(1:4).eq.'-dz=') then ; call read_arg(arg(5:), 0, e_dz, '', 'real', 'dz', eflag)
   else if (arg(1:5).eq.'-min=') then ; call read_arg(arg(6:), 0, e_zmin, '', 'real', 'min', eflag)
   else if (arg(1:5).eq.'-max=') then ; call read_arg(arg(6:), 0, e_zmax, '', 'real', 'max', eflag)
   else ; eflag = .true. ; if (log_errors) write(99,*) "I don't understand the argument: electro "//trim(arg) ; end if

end subroutine read_electro_arg

subroutine read_radial_arg(arg, eflag, log_errors, rad_ws, rad_bins)

   implicit none

   logical(1) :: eflag, log_errors
   character(*) :: arg
   character(5) :: rad_ws(2)
   integer :: rad_bins

   ! if (arg(1:4).eq.'-ws=') then ; call read_arg(arg(5:), rad_ws, 0.0, '', 'int', 'ws', eflag)
   if (arg(1:5).eq.'-ws1=') then ; call read_arg(arg(6:), 0, 0.0, rad_ws(1), 'str', 'ws1', eflag)
   else if (arg(1:5).eq.'-ws2=') then ; call read_arg(arg(6:), 0, 0.0, rad_ws(2), 'str', 'ws2', eflag)
   else if (arg(1:6).eq.'-bins=') then ; call read_arg(arg(7:), rad_bins, 0.0, '', 'int', 'bins', eflag)
   else ; eflag = .true. ; if (log_errors) write(99,*) "I don't understand the argument: radial "//trim(arg) ; end if

end subroutine read_radial_arg

subroutine read_hydration_arg(arg, eflag, log_errors, nh_bins, nh_rcut)

   implicit none

   logical(1) :: eflag, log_errors
   character(*) :: arg
   integer :: nh_bins
   real :: nh_rcut

   if (arg(1:6).eq.'-bins=') then ; call read_arg(arg(7:), nh_bins, 0.0, '', 'int', 'bins', eflag)
   else if (arg(1:5).eq.'-cut=') then ; call read_arg(arg(6:), 0, nh_rcut, '', 'real', 'nhcut', eflag)
   else ; eflag = .true. ; if (log_errors) write(99,*) "I don't understand the argument: hydration "//trim(arg) ; end if

end subroutine read_hydration_arg


subroutine read_rings_input(eflag, r_ns, r_wr, r_ws, r_wh, maxr, maxr_RINGS)

   logical(1) :: in_arg
   character(31) :: buf
   character(5) :: current_arg
   integer :: io, i, j, line_loc

   logical(1) :: eflag
   integer :: maxr, maxr_RINGS, r_ns
   character(*), allocatable :: r_wr(:), r_ws(:)
   integer, allocatable :: r_wh(:,:)


   ! Convert from actual max depth of rings search (MAXR, either 4->6 or 5->7,8,9 up to now...)
   ! To R.I.N.G.S. notation
   if (maxr.eq.6) then ; maxr_RINGS=4 ; else if (maxr.eq.7.or.maxr.eq.8.or.maxr.eq.9) then ; maxr_RINGS=5
   else ; write(99,*) "Values of 6,7,8 or 9 only are allowed for MAXR at the moment..." ; stop ; endif

   ! Read rings input file...
   open(unit=100, file='hin_structure.rings.in', status='old')
   do ; read(100, '(A)', iostat=io) buf ; if (io.ne.0) exit ; r_ns = r_ns+1 ; end do ; close(100)

   if (r_ns.eq.0) then
      write(99,*) "File 'hin_structure.rings.in' not found. Calculating rings for default atoms (SOL:OW, HO4:OW)."
      allocate(r_wr(2), r_ws(2), r_wh(2,4)) ; r_wh = 0 ; r_wr = 'OW'
      r_wr(1) = 'SOL' ; r_wr(2) = 'HO4' ; r_wh(1,1) = 1 ; r_wh(1,2) = 2 ; r_wh(2,1) = 1 ; r_wh(2,2) = 2
   else
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
   end if

end subroutine read_rings_input


subroutine read_b_rcut(eflag, args, b_rcut, npairs)
   implicit none

   logical(1) :: eflag
   integer :: npairs, io, i
   character(*) :: args(npairs)
   real, allocatable :: b_rcut(:)

   do i=1,npairs
      read(args(i), *, iostat=io) b_rcut(i)
      if (io.ne.0) then
         print *, trim(args(i))//" is not a valid value for bonds -rcut!"
         eflag=.true.
      end if
   end do


end subroutine read_b_rcut


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

subroutine set_op_max_cut(switch_qd, switch_qt, q_cut, qd_cut, qt_cut, f_cut, t_rcut, op_max_cut)

   implicit none

   logical(1) :: switch_q(3:6), switch_qt(3:6), switch_qd(3:6)
   real :: q_cut, qd_cut, qt_cut, f_cut, t_rcut, op_max_cut, second_cut=0.0

   if (op_max_cut.lt.q_cut) op_max_cut = q_cut
   if (op_max_cut.lt.f_cut) op_max_cut = f_cut
   if (op_max_cut.lt.t_rcut) op_max_cut = t_rcut
   if (switch_qd(3).or.switch_qd(4).or.switch_qd(6)) second_cut = qd_cut
   if ((second_cut.lt.qt_cut).and.(switch_qd(3).or.switch_qd(4).or.switch_qd(6))) second_cut = qt_cut
   op_max_cut = op_max_cut + second_cut

end subroutine set_op_max_cut

end module MOD_read_input_args

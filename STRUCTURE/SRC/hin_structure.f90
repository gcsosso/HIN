program hin_structure

use, intrinsic :: iso_c_binding, only: C_NULL_CHAR, C_PTR, c_f_pointer
use xtc
use MOD_read_input
use MOD_zdens
use MOD_xyfes
use MOD_rings
use MOD_clusters
use MOD_clathrates
use MOD_bonds
use MOD_electro
use MOD_output
use MOD_order
use MOD_bondorder

implicit none

integer :: NFRAMES, EST_NFRAMES
type(C_PTR) :: OFFSETS

integer, parameter :: dp = kind(1.d0)
integer,parameter :: cart=3, six=6
integer :: NATOMS, STEP, STAT, STAT_OUT, i, j, k, l, m, nz, e_nz, b_bins, nz_bAVE, io, nm
integer ::  ns, idx, nat, dostuff, counter, nl, endf, nxyz, id, ck, ckr, ibin
integer :: per1, per2, per3, per4, per5, per6, kper135, kper246, r13, r15, r24, r26, nper, n_ddc, n_hc
integer :: nsix, r_flag, r_flag2, r_flag3, npairs, npairs_cn, flag, patch, o_nz, f_zbins
integer :: tmplist, ohstride, pmpi, nxy, nsurf, nbulk, nq
integer, allocatable :: n_ws(:), n_r_ws(:), list_ws(:,:), list_r_ws(:,:), r_nper(:), mflag(:), resnum(:)
integer, allocatable :: kto(:), r_color(:), r_array(:), p_rings(:,:,:), C_size(:), C_idx(:,:)
integer :: n_f_ow
integer, allocatable :: list_f_ow(:)
real :: prec, box(cart,cart), box_trans(cart,cart), time, dummyp, lb, ub, icell(cart*cart)
real :: zmin, zmax, dz, rsqdf, posi(cart), posj(cart), xymin, xymax, ddx, ddy, thr
real :: b_zmin, b_zmax, b_dz, b_bmin, b_bmax, rstep, n_ddc_AVE, n_hc_AVE, n_hex_AVE, n_cls_AVE, zop_AVE
real :: n_ddc_AVE_SURF, n_hc_AVE_SURF, n_hex_AVE_SURF, n_ddc_AVE_BULK, n_hc_AVE_BULK, n_hex_AVE_BULK 
real :: ze_AVE, ze_AVE_BULK, ze_AVE_SURF, e_zmin, e_zmax, e_dz, middle, o_dz, hbdist2
real :: f3_imax, f3_cmax, f4_imax, f4_cmin
real :: delta_AVE, delta_AVE_BULK, delta_AVE_SURF, esse_AVE, esse_AVE_BULK, esse_AVE_SURF, rog_AVE, rog_AVE_BULK, rog_AVE_SURF
real, allocatable :: pos(:,:), dens(:,:), zmesh(:), pdbon(:,:,:), stat_nr_AVE(:), xmesh(:), ymesh(:)
real, allocatable :: b_r_cut(:), pdbon_AVE(:,:,:), cn(:,:), cn_AVE(:,:), xydens(:,:,:), stat_nr_HB_AVE(:)
real, allocatable :: d_charge(:), e_zmesh(:), qqq(:), qqq_all(:), mq(:), mq_all(:), w_order(:), o_zmesh(:)
character :: ch
character*3 :: cls_stat
logical(1) :: switch_zdens=.false., switch_cls=.false., switch_bonds=.false., switch_xyfes=.false.
logical(1) :: switch_hw_ex=.true.
logical(1) :: switch_r_idx=.false., switch_electro=.false., switch_f_cls=.false.
character*5, allocatable :: resname(:)
character*4 :: wmol, axis_1, axis_2
character*4, allocatable :: ws(:), sym(:)
character*4, allocatable :: atq(:)
character*100 :: xtcOfile, wformat, natformat, command, plumed_exe
character*100 :: pstring, pstring_C, command1, command2, fcommand, vmd_exe, buffer
type(C_PTR) :: xd_c, xd_c_out
type(xdrfile), pointer :: xd, xd_out
logical(1) :: ex, proc, cknn
real(dp), allocatable :: qlb_io(:)

integer, parameter :: ARG_LEN=127

! TRAJECTORY
character(ARG_LEN) :: sfile='md.gro', tfile='traj.xtc'
integer :: fframe=0, lframe=-1, stride=1
logical(1) :: switch_outxtc=.true., switch_progress=.false.

! ORDER
logical(1) :: switch_op=.false., switch_q(3:6)=.false., switch_qd(3:6)=.false., switch_qt(3:6)=.false.
logical(1) :: switch_t4=.false., switch_f(3:4)=.false., switch_th=.false.
character(7) :: filter='none'
character(5) :: op_species='OW'
character(3) :: switch_water='mol'
real :: filt_min=0.0, filt_max=1.0, q_cut=0.35, qd_cut=0.35, qt_cut=0.35, f_cut=0.35
integer :: max_shell=30

! RINGS
logical(1) :: switch_rings=.false., switch_r_split=.false., switch_hbck=.false., switch_hex=.false.
logical(1) :: switch_r_cls=.false., switch_cages=.false., switch_ffss=.false.
character(ARG_LEN) :: rings_exe=''
character(3) :: r_cls_W=''
real :: r_split=0.0, r_cut=0.32, hbdist=0.22, hbangle=20.0, a_thr=30.0, thrS=1.5, thrSS=1.7
integer :: maxr=9, maxr_RINGS, wcol=6, r_ns=0
character(5), allocatable :: r_wr(:), r_ws(:)
integer, allocatable :: r_wh(:,:)

!REPLACE THIS!
ns = 1
allocate(n_ws(ns), ws(ns))
ws(1) = 'OW'

! Open the .log file
open(unit=99, file='hin_structure.log', status='unknown')

call read_input(ARG_LEN,sfile,tfile,fframe,lframe,stride,switch_outxtc, switch_progress, &
                switch_op,switch_q,switch_qd,switch_qt,switch_t4,switch_f,switch_th, &
                filter,filt_min,filt_max,q_cut,qd_cut,qt_cut,f_cut,max_shell,op_species, &
                switch_rings, switch_r_split, switch_hbck, switch_hex, switch_r_cls, switch_cages, switch_ffss, &
                rings_exe, r_cls_W, r_split, r_cut, hbdist, hbangle, a_thr, thrS, thrSS, maxr, maxr_RINGS, wcol, &
                r_ns, r_wr, r_ws, r_wh)

if (lframe.eq.-1) then
   STAT=read_xtc_n_frames(trim(adjustl(tfile))//C_NULL_CHAR, NFRAMES, EST_NFRAMES, OFFSETS)
   if (STAT.ne.0) then ; write(99,*) "Something is wrong with the trajectory file..." ; stop
   else ; lframe = NFRAMES-1 ; end if
end if

hbdist2 = hbdist**2.0

call read_gro(sfile,nat,sym,list_ws,list_r_ws,r_color,kto,n_ws,switch_hw_ex,switch_rings,r_ns,r_ws,r_wr,n_r_ws, &
              natformat,ns,resnum,resname,idx,dummyp,ws,list_f_ow,n_f_ow,switch_op)

allocate(qlb_io(n_f_ow))

!! JPCL stuff : read the flags that tell you whether a conf. is surviving or dying
!open(unit=877, file='flags.dat', status='old')
!nm=0
!do
!  read(877,*,iostat=io) buffer
!  if (io/=0) exit
!  nm=nm+1
!end do
!rewind(877)
!allocate(mflag(nm))
!do i=1,nm
!   read(877,*) mflag(i)
!end do
!close(877)

! Set the averages to zero !
nz_bAVE=0; n_ddc_AVE=0.0; n_hc_AVE=0.0; n_hex_AVE=0.0; n_cls_AVE=0.0; zop_AVE=0.0; 
n_ddc_AVE_SURF=0.0; n_hc_AVE_SURF=0.0; n_hex_AVE_SURF=0.0; n_ddc_AVE_BULK=0.0; 
n_hc_AVE_BULK=0.0; n_hex_AVE_BULK=0.0; ze_AVE=0.0; ze_AVE_BULK=0.0;  ze_AVE_SURF=0.0; 
delta_AVE=0.0;  delta_AVE_BULK=0.0;  delta_AVE_SURF=0.0; esse_AVE=0.0; esse_AVE_BULK=0.0; 
esse_AVE_SURF=0.0; rog_AVE=0.0; rog_AVE_BULK=0.0; rog_AVE_SURF=0.0; 

if (switch_zdens) call zdens_alloc(nz,zmax,zmin,dz,dens,zmesh,ns)

call read_first_xtc(tfile,switch_outxtc,xtcOfile,STAT,NATOMS,nat,xd_c,xd,xd_c_out,xd_out,STEP,time,box_trans,pos,prec,icell,cart)

if (switch_rings.and.switch_r_idx) call read_cls_idx(lframe,fframe,stride,C_size,C_idx,nat)

if (switch_xyfes) call xyfes_alloc(nxy,xymax,xymin,ddx,ddy,xydens,xmesh,ymesh,ns,icell,cart)

if (switch_rings) then
        call rings_alloc(switch_rings,switch_cages,switch_hex,switch_outxtc, &
                 stat_nr_AVE,maxr,n_ddc_AVE,n_hc_AVE,n_hex_AVE,switch_r_idx, &
                 switch_r_cls,r_cls_W,nsurf,nbulk,n_ddc_AVE_SURF,n_hc_AVE_SURF,n_hex_AVE_SURF, &
                 n_ddc_AVE_BULK,n_hc_AVE_BULK,n_hex_AVE_BULK,delta_AVE,delta_AVE_BULK, &
                 delta_AVE_SURF,esse_AVE,esse_AVE_BULK,esse_AVE_SURF,rog_AVE,rog_AVE_BULK,rog_AVE_SURF, &
                 ze_AVE,ze_AVE_BULK,ze_AVE_SURF,stat_nr_HB_AVE,switch_hbck)
end if

if (switch_cls) call clusters_alloc(switch_outxtc,ns,ws,switch_hw_ex,ohstride,n_ws,list_ws,sfile,vmd_exe,n_cls_AVE)

if (switch_f(3).or.switch_f(4)) call clathrates_alloc(switch_f,switch_f_cls)

if (switch_bonds) call bonds_alloc(nz_bAVE,b_zmax,b_zmin,b_dz,pdbon_AVE,npairs,b_bins,cn_AVE,npairs_cn)

if (switch_electro.or.switch_th) then ! we need masses for order as well - COM removal...
   call electro_alloc(e_zmin,e_zmax,e_dz,d_charge, &
        e_zmesh,e_nz,atq,qqq,nq,qqq_all,nat,sym,middle,pos,cart,mq,mq_all)
end if

!if (trim(adjustl(switch_th)).eq.'yes'.and.trim(adjustl(switch_water)).eq.'yes') then
if (switch_th) call order_alloc(o_nz,filt_max,filt_min,o_dz,w_order,o_zmesh,switch_water)

if (switch_t4) call bondorder_t4_alloc()
if (switch_q(3).or.switch_qd(3).or.switch_qt(3)) call bondorder_alloc(3)
if (switch_q(4).or.switch_qd(4).or.switch_qt(4)) call bondorder_alloc(4)
if (switch_q(6).or.switch_qd(6).or.switch_qt(6)) call bondorder_alloc(6)

! Read the whole thing
counter=0
dostuff=0
call progress(0.0, switch_progress)

do while ( STAT==0 )
   if (mod(counter,stride).eq.0.and.counter.ge.fframe.and.counter.le.lframe) then
      write(99,'(a,f18.6,a,i0,a,i0)') " Time (ps): ", time, "  Step: ", STEP, " Frame: ", counter
      dostuff=dostuff+1

      ! Write .xtc...
      if (switch_outxtc) STAT_OUT=write_xtc(xd_out,NATOMS,STEP,time,box_trans,pos,prec)

      ! Number density profile along z...
      if (switch_zdens) call zdens(ns,n_ws,nz,zmin,dz,pos,cart,list_ws,dens)

      ! 2D FES in the xy plane...
      if (switch_xyfes) call xyfes(ns,n_ws,nxy,xymin,ddx,ddy,pos,cart,list_ws,xydens,xymax,icell)

      ! Rings statistics...  
      if (switch_rings) then
          call rings(kto,r_ns,r_wh,n_r_ws,pos,cart,list_r_ws,filt_min,filt_max,sym,resname,rings_exe,r_color,time,STEP, &
                     counter,natformat,nat,icell,r_cut,n_ddc_AVE,n_hc_AVE,a_thr,maxr,maxr_RINGS,switch_r_split,r_split, &
                     switch_cages,stat_nr_AVE,switch_hex,n_hex_AVE,wcol,box_trans,switch_r_cls,r_cls_W, &
                     patch,switch_r_idx,C_size,C_idx,switch_ffss,thrS,nsurf,nbulk,n_ddc_AVE_SURF, &
                     n_hc_AVE_SURF,n_hex_AVE_SURF, &
                     n_ddc_AVE_BULK,n_hc_AVE_BULK,n_hex_AVE_BULK, &
                     delta_AVE,delta_AVE_BULK,delta_AVE_SURF,esse_AVE,esse_AVE_BULK,esse_AVE_SURF, &
                     rog_AVE,rog_AVE_BULK,rog_AVE_SURF,ze_AVE,ze_AVE_BULK,ze_AVE_SURF,mflag,switch_hbck, &
                     hbdist,hbdist2,hbangle,stat_nr_HB_AVE,thrSS)

      end if

      ! Ice-like clusters...
      if (switch_cls) call clusters(plumed_exe,ns,ws,list_ws,NATOMS,STEP,time, &
                                    box_trans,pos,prec,cart,pmpi,cls_stat,r_color,natformat,nat,n_cls_AVE)
      
      ! Clathrates...
      if (switch_f(3).or.switch_f(4)) then
        call clathrates(switch_f,filt_min,filt_max,f_cut,n_f_ow,list_f_ow,counter, &
                        time,cart,icell,pos,nat,natformat,f_zbins,switch_f_cls,f3_imax,f3_cmax,f4_imax,f4_cmin)
      end if
      
      ! Bonds statistics...
      if (switch_bonds) call bonds(b_zmin,b_zmax,b_dz,pos,icell,pdbon,cart,ns,n_ws,list_ws, &
                                   ws,b_r_cut,npairs,b_bins,b_bmin,b_bmax,cn,npairs_cn,cn_AVE,pdbon_AVE)


      ! Electrostatic...
      if (switch_electro) call electro(e_zmin,e_zmax,e_dz,nat,e_nz, &
                                       pos,d_charge,cart,sym,atq,qqq,nq,qqq_all,middle,mq,mq_all)

      ! Ordering...
      if (switch_th) then
         call order(o_nz,filt_max,filt_min,o_dz,w_order,o_zmesh,nat,pos, &
              mq_all,cart,middle,switch_water,sym,wmol,resname, &
              resnum,axis_1,axis_2,zop_AVE,natformat,icell)
      end if
		
      ! Q Ordering...
      if (switch_q(3).or.switch_qd(3).or.switch_qt(3)) then
          call bondorder(3,filt_min,filt_max,q_cut,qd_cut,qt_cut,counter,list_f_ow,n_f_ow,max_shell, &
                         time,cart,icell,pos,nat,natformat,sym,switch_q(3),switch_qd(3),switch_qt(3),switch_t4,qlb_io)
      end if
      if (switch_q(4).or.switch_qd(4).or.switch_qt(4)) then
          call bondorder(4,filt_min,filt_max,q_cut,qd_cut,qt_cut,counter,list_f_ow,n_f_ow,max_shell, &
                         time,cart,icell,pos,nat,natformat,sym,switch_q(4),switch_qd(4),switch_qt(4),switch_t4,qlb_io)
      end if
      if (switch_q(6).or.switch_qd(6).or.switch_qt(6)) then
          call bondorder(6,filt_min,filt_max,q_cut,qd_cut,qt_cut,counter,list_f_ow,n_f_ow,max_shell, &
                         time,cart,icell,pos,nat,natformat,sym,switch_q(6),switch_qd(6),switch_qt(6),switch_t4,qlb_io)
      end if
      
   end if
   
   counter=counter+1
   call progress(real(counter-fframe)/real(lframe-fframe+1), switch_progress)
   if (counter.gt.lframe) exit
   ! Read .xtc frame...
   STAT=read_xtc(xd,NATOMS,STEP,time,box_trans,pos,prec)
end do

! Output...
call output(dostuff,lframe,fframe,stride,switch_outxtc,ns,ws,n_ws,zmesh,dens,nz,dz,box_trans, &
            switch_rings,r_ns,r_ws,n_r_ws,maxr,stat_nr_AVE,switch_cages,n_ddc_AVE,n_hc_AVE, &
            switch_hex,n_hex_AVE,switch_bonds,npairs,nz_bAVE,b_zmin,b_dz,b_bins,b_bmax,b_bmin, &
            pdbon_AVE,npairs_cn,cn_AVE,switch_cls,n_cls_AVE,cart,switch_zdens, &
            switch_r_cls,r_cls_W,switch_xyfes,xydens,xymax,xymin,nxy,xmesh,ymesh,nsurf,nbulk,n_ddc_AVE_SURF, &
            n_hc_AVE_SURF,n_hex_AVE_SURF,n_ddc_AVE_BULK,n_hc_AVE_BULK,n_hex_AVE_BULK,switch_ffss, &
            delta_AVE,delta_AVE_BULK,delta_AVE_SURF,esse_AVE,esse_AVE_BULK,esse_AVE_SURF, &
            rog_AVE,rog_AVE_BULK,rog_AVE_SURF,ze_AVE,ze_AVE_BULK,ze_AVE_SURF,d_charge, &
            switch_electro,e_nz,e_zmesh,switch_th,switch_water,o_nz,o_zmesh,w_order,zop_AVE,stat_nr_HB_AVE,switch_hbck)

STAT=xdrfile_close(xd)

! At a certain point I shuld consider to deallocate this stuff...
!deallocate(pos,sym,ws,n_ws,dens)

! Cleaning up ! If you need some stuff for debugging purposes, just comment the following line
if (switch_rings.and.switch_r_idx) call system('mv idx.dat idx.tempo')
call system('rm -r -f rings.in tmp* rstat data bonds Walltime options log conf.pdb vmd.log')
call system('rm -r -f rings.out plumed.xtc r3-5.dat r4-5.dat r5-5.dat r6-5.dat')
if (switch_rings.and.switch_r_idx) call system('mv idx.tempo idx.dat')

end program hin_structure

subroutine progress(j)

   implicit none
   
   logical(1) :: switch_progress
   integer :: k
   real :: j
   character(57) :: bar="???% |                                                  |"
   
   if (switch_progress) then
      write(unit=bar(1:3),fmt="(i3)") int(100*j)
      do k=1, int(50*j)
         bar(6+k:6+k)="*"
      end do
      ! print the progress bar.
      write(unit=6,fmt="(a1,a60,$)") char(13), bar
   end if
end subroutine progress
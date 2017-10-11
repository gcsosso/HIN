program hin_structure

use, intrinsic :: iso_c_binding, only: C_NULL_CHAR, C_PTR, c_f_pointer
use xtc
use MOD_read_input
use MOD_zdens
use MOD_xyfes
use MOD_rings
use MOD_clusters
use MOD_bonds
use MOD_electro
use MOD_output
use MOD_order

implicit none

integer,parameter :: cart=3, six=6
integer :: NATOMS, STEP, STAT, STAT_OUT, i, j, k, l, m, fframe, stride, lframe, nz, e_nz, b_bins, nz_bAVE, io, nm
integer ::  eflag, ns, r_ns, idx, nat, dostuff, counter, nl, endf, nxyz, id, ck, ckr, ibin
integer :: per1, per2, per3, per4, per5, per6, kper135, kper246, r13, r15, r24, r26, nper, n_ddc, n_hc
integer :: nsix, r_flag, r_flag2, r_flag3, npairs, npairs_cn, flag, patch, o_nz
integer :: maxr, maxr_RINGS, wcol, tmplist, ohstride, pmpi, nxy, nsurf, nbulk, nq
integer, allocatable :: n_ws(:), n_r_ws(:), list_ws(:,:), list_r_ws(:,:), r_nper(:), mflag(:), resnum(:)
integer, allocatable :: kto(:), r_color(:), r_array(:), p_rings(:,:,:), C_size(:), C_idx(:,:)
real :: prec, box(cart,cart), box_trans(cart,cart), time, dummyp, lb, ub, icell(cart*cart)
real :: zmin, zmax, r_zmin, r_zmax, dz, rcut, rsqdf, posi(cart), posj(cart), xymin, xymax, ddx, ddy, thr, thrS, thrSS
real :: b_zmin, b_zmax, b_dz, b_bmin, b_bmax, rstep, a_thr, n_ddc_AVE, n_hc_AVE, n_hex_AVE, n_cls_AVE, zop_AVE
real :: n_ddc_AVE_SURF, n_hc_AVE_SURF, n_hex_AVE_SURF, n_ddc_AVE_BULK, n_hc_AVE_BULK, n_hex_AVE_BULK 
real :: ze_AVE, ze_AVE_BULK, ze_AVE_SURF, e_zmin, e_zmax, e_dz, middle, o_zmax, o_zmin, o_dz, hbdist, hbangle
real :: delta_AVE, delta_AVE_BULK, delta_AVE_SURF, esse_AVE, esse_AVE_BULK, esse_AVE_SURF, rog_AVE, rog_AVE_BULK, rog_AVE_SURF
real, allocatable :: pos(:,:), dens(:,:), zmesh(:), pdbon(:,:,:), stat_nr_AVE(:), xmesh(:), ymesh(:)
real, allocatable :: b_rcut(:), pdbon_AVE(:,:,:), cn(:,:), cn_AVE(:,:), xydens(:,:,:), stat_nr_HB_AVE(:)
real, allocatable :: d_charge(:), e_zmesh(:), qqq(:), qqq_all(:), mq(:), mq_all(:), w_order(:), o_zmesh(:)
character :: ch
character*3 :: outxtc, hw_ex, switch_zdens, switch_rings, switch_cls, switch_bonds, switch_xyfes
character*3 :: switch_hex, switch_r_cls, r_cls_W, switch_cages, cls_stat, switch_r_idx, switch_ffss
character*3 :: switch_electro, switch_order, switch_water, switch_hbck
character*5, allocatable :: resname(:)
character*4 :: wmol, axis_1, axis_2
character*4, allocatable :: ws(:), r_ws(:), sym(:)
character*4, allocatable :: atq(:)
character*100 :: sfile, tfile, xtcOfile, wformat, natformat, rings_exe, command, plumed_exe
character*100 :: pstring, pstring_C, command1, command2, fcommand, vmd_exe, buffer
type(C_PTR) :: xd_c, xd_c_out
type(xdrfile), pointer :: xd, xd_out
logical :: ex, proc, cknn

! Open the .log file
open(unit=99, file='hin_structure.log', status='unknown')

call read_input(eflag,sfile,tfile,fframe,stride,lframe,outxtc,hw_ex,switch_zdens,ns,ws,n_ws,zmin,zmax,dz, & 
                switch_rings,rings_exe,r_zmin,r_zmax,r_ns,r_ws,n_r_ws,rcut,switch_cls,plumed_exe, &
                switch_bonds,b_zmin,b_zmax,b_dz,b_rcut,npairs,b_bins,b_bmin,b_bmax,npairs_cn,maxr, &
                switch_hex,switch_r_cls,r_cls_W,a_thr,maxr_RINGS,switch_cages,wcol,ohstride, &
                vmd_exe,pmpi,cls_stat,switch_xyfes,xymin,xymax,nxy,switch_r_idx,switch_ffss,thrS, &
                switch_electro,e_zmin,e_zmax,e_dz,switch_order,wmol,axis_1,axis_2,o_zmin, &
                o_zmax,o_dz,switch_water,switch_hbck,hbdist,hbangle,thrSS)

call read_gro(sfile,nat,sym,list_ws,list_r_ws,r_color,kto,n_ws,hw_ex,switch_rings,r_ns,r_ws,n_r_ws, &
              natformat,ns,resnum,resname,idx,dummyp,ws)


!! JPCL stuff : read the flags that tell you whether a conf. is surviving or dying
!open(unit=877, file='flags.dat', status='old')
!nm=0
!do
!  read(877,*,iostat=io) buffer
!  if (io/=0) exit
!  nm=nm+1
!enddo
!rewind(877)
!allocate(mflag(nm))
!do i=1,nm
!   read(877,*) mflag(i)
!enddo
!close(877)


if (trim(adjustl(switch_zdens)).eq.'yes') then
   call zdens_alloc(nz,zmax,zmin,dz,dens,zmesh,ns)
endif

call read_first_xtc(tfile,outxtc,xtcOfile,STAT,NATOMS,nat,xd_c,xd,xd_c_out,xd_out,STEP,time,box_trans,pos,prec,icell,cart)

if (trim(adjustl(switch_rings)).eq.'yes'.and.trim(adjustl(switch_r_idx)).eq.'yes') then
   call read_cls_idx(lframe,fframe,stride,C_size,C_idx,nat)
endif

if (trim(adjustl(switch_xyfes)).eq.'yes') then
   call xyfes_alloc(nxy,xymax,xymin,ddx,ddy,xydens,xmesh,ymesh,ns,icell,cart)
endif

if (trim(adjustl(switch_rings)).eq.'yes') then
        call rings_alloc(switch_rings,switch_cages,switch_hex,outxtc, &
                 stat_nr_AVE,maxr,n_ddc_AVE,n_hc_AVE,n_hex_AVE, &
                 switch_r_cls,r_cls_W,nsurf,nbulk,n_ddc_AVE_SURF,n_hc_AVE_SURF,n_hex_AVE_SURF, &
                 n_ddc_AVE_BULK,n_hc_AVE_BULK,n_hex_AVE_BULK,delta_AVE,delta_AVE_BULK, &
                 delta_AVE_SURF,esse_AVE,esse_AVE_BULK,esse_AVE_SURF,rog_AVE,rog_AVE_BULK,rog_AVE_SURF, &
                 ze_AVE,ze_AVE_BULK,ze_AVE_SURF,stat_nr_HB_AVE)
endif

if (trim(adjustl(switch_cls)).eq.'yes') then
   call clusters_alloc(switch_cls,outxtc,ns,ws,hw_ex,ohstride,n_ws,list_ws,sfile,vmd_exe,n_cls_AVE)
   endif

if (trim(adjustl(switch_bonds)).eq.'yes') then
   call bonds_alloc(nz_bAVE,b_zmax,b_zmin,b_dz,pdbon_AVE,npairs,b_bins,cn_AVE,npairs_cn)
endif

if (trim(adjustl(switch_electro)).eq.'yes'.or.(trim(adjustl(switch_order)).eq.'yes')) then ! we need masses for order as well - COM removal...
   call electro_alloc(e_zmin,e_zmax,e_dz,d_charge, &
        e_zmesh,e_nz,atq,qqq, &
        nq,qqq_all,nat,sym,middle,pos,cart,mq,mq_all)
endif

if (trim(adjustl(switch_order)).eq.'yes'.and.trim(adjustl(switch_water)).eq.'yes') then
   call order_alloc(o_nz,o_zmax,o_zmin,o_dz,w_order,o_zmesh)
endif


! Read the whole thing
counter=0
dostuff=0

do while ( STAT==0 )
   if (mod(counter,stride).eq.0.and.counter.ge.fframe.and.counter.le.lframe) then
      write(99,'(a,f18.6,a,i0,a,i0)') " Time (ps): ", time, "  Step: ", STEP, " Frame: ", counter
      dostuff=dostuff+1

      ! Write .xtc...
      if (trim(adjustl(outxtc)).eq.'yes') then
         STAT_OUT=write_xtc(xd_out,NATOMS,STEP,time,box_trans,pos,prec) 
      endif

      ! Number density profile along z...
      if (trim(adjustl(switch_zdens)).eq.'yes') then
         call zdens(switch_zdens,ns,n_ws,nz,zmin,dz,pos,cart,list_ws,dens)
      endif

      ! 2D FES in the xy plane...
      if (trim(adjustl(switch_xyfes)).eq.'yes') then
         call xyfes(switch_xyfes,ns,n_ws,nxy,xymin,ddx,ddy,pos,cart,list_ws,xydens,xymax,icell)
      endif

      ! Rings statistics...  
      if (trim(adjustl(switch_rings)).eq.'yes') then
          call rings(kto,r_ns,n_r_ws,pos,cart,list_r_ws,r_zmin,r_zmax,sym,rings_exe,r_color,time,STEP, &
                     counter,natformat,nat,icell,rcut,n_ddc_AVE,n_hc_AVE,a_thr,maxr,maxr_RINGS, &
                     switch_cages,stat_nr_AVE,switch_hex,n_hex_AVE,wcol,box_trans,switch_r_cls,r_cls_W, &
                     patch,switch_r_idx,C_size,C_idx,switch_ffss,thrS,nsurf,nbulk,n_ddc_AVE_SURF, &
                     n_hc_AVE_SURF,n_hex_AVE_SURF, &
                     n_ddc_AVE_BULK,n_hc_AVE_BULK,n_hex_AVE_BULK, &
                     delta_AVE,delta_AVE_BULK,delta_AVE_SURF,esse_AVE,esse_AVE_BULK,esse_AVE_SURF, &
                     rog_AVE,rog_AVE_BULK,rog_AVE_SURF,ze_AVE,ze_AVE_BULK,ze_AVE_SURF,mflag,switch_hbck, &
                     hbdist,hbangle,stat_nr_HB_AVE,thrSS)

      endif

      ! Ice-like clusters...
      if (trim(adjustl(switch_cls)).eq.'yes') then 
         call clusters(outxtc,plumed_exe,ns,ws,list_ws,NATOMS,STEP,time, &
                       box_trans,pos,prec,cart,pmpi,cls_stat,r_color,natformat,nat,n_cls_AVE)
      endif
      
      ! Bonds statistics...
      if (trim(adjustl(switch_bonds)).eq.'yes') then 
         call bonds(b_zmin,b_zmax,b_dz,pos,icell,pdbon,cart,ns,n_ws,list_ws, &
                    ws,b_rcut,npairs,b_bins,b_bmin,b_bmax,cn,npairs_cn,cn_AVE,pdbon_AVE)
      endif


      ! Electrostatic...
      if (trim(adjustl(switch_electro)).eq.'yes') then
         call electro(e_zmin,e_zmax,e_dz,nat,e_nz, &
              pos,d_charge,cart,sym,atq,qqq,nq,qqq_all,middle,mq,mq_all)
      endif

      ! Ordering...
      if (trim(adjustl(switch_order)).eq.'yes') then
         call order(o_nz,o_zmax,o_zmin,o_dz,w_order,o_zmesh,nat,pos, &
              mq_all,cart,middle,switch_water,sym,wmol,resname, &
              resnum,axis_1,axis_2,zop_AVE)
      endif
      
   endif
   counter=counter+1
   if (counter.gt.lframe) exit
   ! Read .xtc frame...
   STAT=read_xtc(xd,NATOMS,STEP,time,box_trans,pos,prec)
enddo

! Output...
call output(dostuff,lframe,fframe,stride,outxtc,ns,ws,n_ws,zmesh,dens,nz,dz,box_trans, &
            switch_rings,r_ns,r_ws,n_r_ws,maxr,stat_nr_AVE,switch_cages,n_ddc_AVE,n_hc_AVE, &
            switch_hex,n_hex_AVE,switch_bonds,npairs,nz_bAVE,b_zmin,b_dz,b_bins,b_bmax,b_bmin, &
            pdbon_AVE,npairs_cn,cn_AVE,switch_cls,n_cls_AVE,cart,switch_zdens, &
            switch_r_cls,r_cls_W,switch_xyfes,xydens,xymax,xymin,nxy,xmesh,ymesh,nsurf,nbulk,n_ddc_AVE_SURF, &
            n_hc_AVE_SURF,n_hex_AVE_SURF,n_ddc_AVE_BULK,n_hc_AVE_BULK,n_hex_AVE_BULK,switch_ffss, &
            delta_AVE,delta_AVE_BULK,delta_AVE_SURF,esse_AVE,esse_AVE_BULK,esse_AVE_SURF, &
            rog_AVE,rog_AVE_BULK,rog_AVE_SURF,ze_AVE,ze_AVE_BULK,ze_AVE_SURF,d_charge, &
            switch_electro,e_nz,e_zmesh,switch_order,switch_water,o_nz,o_zmesh,w_order,zop_AVE,stat_nr_HB_AVE,switch_hbck)

STAT=xdrfile_close(xd)

! At a certain point I shuld consider to deallocate this stuff...
!deallocate(pos,sym,ws,n_ws,dens)

! Cleaning up ! If you need some stuff for debugging purposes, just comment the following line
if (trim(adjustl(switch_rings)).eq.'yes'.and.trim(adjustl(switch_r_idx)).eq.'yes') then
   command='mv idx.dat idx.tempo'
   call system(command)
endif
command='rm -r -f rings.in tmp* rstat data bonds Walltime options log conf.pdb vmd.log rings.out plumed.xtc r*-5.dat'
call system(command)
if (trim(adjustl(switch_rings)).eq.'yes'.and.trim(adjustl(switch_r_idx)).eq.'yes') then
   command='mv idx.tempo idx.dat'
   call system(command)
endif

end program hin_structure

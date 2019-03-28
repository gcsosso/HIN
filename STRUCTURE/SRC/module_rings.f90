module MOD_rings

contains

subroutine rings_alloc(switch_rings,switch_cages,switch_hex,outxtc, &
                       stat_nr_AVE,maxr,n_ddc_AVE,n_hc_AVE,n_hex_AVE, &
                       switch_r_cls,r_cls_W,nsurf,nbulk,n_ddc_AVE_SURF,n_hc_AVE_SURF,n_hex_AVE_SURF, &
                       n_ddc_AVE_BULK,n_hc_AVE_BULK,n_hex_AVE_BULK, &
                       delta_AVE,delta_AVE_BULK,delta_AVE_SURF,esse_AVE,esse_AVE_BULK,esse_AVE_SURF, &
                       rog_AVE,rog_AVE_BULK,rog_AVE_SURF,ze_AVE,ze_AVE_BULK,ze_AVE_SURF,stat_nr_HB_AVE,switch_hbck)

implicit none

! Local
character*100 :: command

! Arguments
integer :: maxr, nsurf, nbulk
real :: n_ddc_AVE, n_hc_AVE, n_hex_AVE, n_ddc_AVE_SURF, n_hc_AVE_SURF, n_hex_AVE_SURF
real :: n_ddc_AVE_BULK, n_hc_AVE_BULK, n_hex_AVE_BULK
real :: ze_AVE,ze_AVE_BULK,ze_AVE_SURF
real :: delta_AVE, delta_AVE_BULK, delta_AVE_SURF, esse_AVE, esse_AVE_BULK, esse_AVE_SURF, rog_AVE, rog_AVE_BULK, rog_AVE_SURF
real, allocatable :: stat_nr_AVE(:), stat_nr_HB_AVE(:)
character*3 :: switch_rings, switch_cages, switch_hex, outxtc, switch_r_cls, r_cls_W, switch_hbck

! If we're doing rings, make the tmp dir and open the output files...
if (trim(adjustl(switch_rings)).eq.'yes') then
   command="rm -r -f data ; mkdir data"
   call system(command)
   if (trim(adjustl(switch_cages)).eq.'yes') then
      write(99,*) "We have looked into DDCs and HCs as well. See hin_structure.out.rings.cages"
      open(unit=103, file='hin_structure.out.rings.cages', status='unknown')
      write(103,*) "# Time [ps] | N. of 6-membered rings | N. of DDC cages | N. of HC cages"
   endif
   if (trim(adjustl(switch_hex)).eq.'yes') then
      write(99,*) "We have looked into hexagonal rings as well. See hin_structure.out.rings.hex"
      open(unit=108, file='hin_structure.out.rings.hex', status='unknown')
      write(108,*) "! Time [ps] | N. of 6-membered rings | N. of proper hexagonal rings"
   endif
   close(108)
   if (trim(adjustl(outxtc)).eq.'no') then
      write(99,*) "hin_structure.out.rings.color could not match the .xtc traj!"
      write(99,*) "switch OUTXTC to yes!"
      stop
   endif
   open(unit=104, file='hin_structure.out.rings.color', status='unknown')
   open(unit=107, file='hin_structure.out.rings.stats', status='unknown')
   write(107,*) "# Time [ps] | N. of n-membered rings (from 3 to n (max=9))"
   if (trim(adjustl(switch_hbck)).eq.'yes') then
      open(unit=307, file='hin_structure.out.rings.stats.HB', status='unknown')
      write(307,*) "# Time [ps] | N. of ** HB ** n-membered rings (from 3 to n (max=9))"
   endif
   ! Cluster hexagonal rings, e.g. to find the largest patch of hexagonal rings sitting on top of the surface
   if (trim(adjustl(switch_r_cls)).eq.'yes') then
      if (trim(adjustl(r_cls_W)).eq.'CLA') then
         open(unit=210, file='hin_structure.out.rings.clath', status='unknown')
         open(unit=211, file='hin_structure.out.rings.clath.color', status='unknown')
         write(210,*) "# Time [ps] | N. of 555 partcages | N. of 655 partcages"
      else if (trim(adjustl(r_cls_W)).ne.'SIX') then
         write(99,*) "Sorry mate, I can do only six membered rings at the moment..."
      else
      open(unit=207, file='hin_structure.out.rings.stats.patch', status='unknown')
      open(unit=208, file='hin_structure.out.rings.color.patch', status='unknown')
      endif
   endif 
   ! Allocate the average...
   allocate(stat_nr_AVE(maxr-2),stat_nr_HB_AVE(maxr-2))
   stat_nr_AVE(:)=0.0
   stat_nr_HB_AVE(:)=0.0
   n_ddc_AVE=0.0
   n_hc_AVE=0.0
   n_hex_AVE=0.0
   nsurf=0 
   nbulk=0
   n_ddc_AVE_BULK=0.0; n_hc_AVE_BULK=0.0; n_hex_AVE_BULK=0.0
   n_ddc_AVE_SURF=0.0; n_hc_AVE_SURF=0.0; n_hex_AVE_SURF=0.0
   delta_AVE=0.0; delta_AVE_BULK=0.0; delta_AVE_SURF=0.0
   esse_AVE=0.0; esse_AVE_BULK=0.0; esse_AVE_SURF=0.0
   rog_AVE=0.0; rog_AVE_BULK=0.0; rog_AVE_SURF=0.0
   ze_AVE=0.0; ze_AVE_BULK=0.0; ze_AVE_SURF=0.0
   
endif

end subroutine rings_alloc

subroutine rings(kto,r_ns,n_r_ws,pos,cart,list_r_ws,r_zmin,r_zmax, &
                 sym,rings_exe,r_color,time,STEP,counter,natformat, &
                 nat,icell,rcut,n_ddc_AVE,n_HC_AVE,a_thr,maxr,maxr_RINGS, &
                 switch_cages,stat_nr_AVE,switch_hex,n_hex_AVE,wcol,box_trans,switch_r_cls,r_cls_W, &
                 patch,switch_r_idx,C_size,C_idx,switch_ffss,thrS,nsurf,nbulk,n_ddc_AVE_SURF, &
                 n_hc_AVE_SURF,n_hex_AVE_SURF,n_ddc_AVE_BULK,n_hc_AVE_BULK,n_hex_AVE_BULK, &
                 delta_AVE,delta_AVE_BULK,delta_AVE_SURF,esse_AVE,esse_AVE_BULK,esse_AVE_SURF, &
                 rog_AVE,rog_AVE_BULK,rog_AVE_SURF,ze_AVE,ze_AVE_BULK,ze_AVE_SURF,mflag,switch_hbck, &
                 hbdist,hbangle,stat_nr_HB_AVE,thrSS)

use dfs
use MOD_vector3

implicit none

! Local
integer, parameter :: osix=6
integer :: i, j, k, l, m, nxyz, nl, endf, iostat, id, nsix, n_hex, cart
integer :: r_flag, r_flag2, r_flag3, tr(osix), tr6(osix), nper, ckr, ck
integer :: per1, per2, per3, per4, per5, per6, kper135, kper246, r13
integer :: r15, r24, r26, n_ddc, n_hc, maxr, maxr_RINGS, info
integer :: nmem, ti, tj, tk, tri, kr, surfF, lwork, dummy
! JPCL
integer :: ddc_bulk_madeit, hc_bulk_madeit, ddc_bulk_dead, hc_bulk_dead, ddc_surf_madeit, hc_surf_madeit, ddc_surf_dead, hc_surf_dead
! JPCL
integer, allocatable :: w_rings(:,:),r_array(:),p_rings(:,:,:),r_nper(:)
integer, allocatable :: stat_nr(:), stat_nr_HB(:)
real :: posi(cart), posj(cart), posk(cart), xdf, ydf, zdf, dist, a_thr
real :: pol1(cart), pol2(cart), rcm(cart), tin(cart,cart), v1(cart), v2(cart), v1m, v2m
real :: rij(cart), rkj(cart), rij_M , rkj_M, hex_angle, z_ext
real, parameter :: rad2deg=57.2958, ref_angle=120.0, pi=4.0*atan(1.0)
double precision, allocatable :: work(:)
double precision :: mtemp(cart,cart), eigen(cart), delta, esse, rog, trt, trt2, lambda, delta1
character*100 :: command, command2, rst, rst2, fcommand, stat_format, arname
logical :: cknn, exist
type(ragged_array) :: stat_wr

! Arguments
integer :: r_ns, STEP, counter, six, nat, wcol, nsurf, nbulk, hbflag(cart*cart) ! no more than 9-membered rings in any case...
integer, allocatable :: kto(:), n_r_ws(:), list_r_ws(:,:), r_color(:), C_size(:), C_idx(:,:), mflag(:)
real :: r_zmin, r_zmax, time, icell(cart*cart), rcut, thrS, minz, cell_ORTO(cart), maxz, thrSS
real ::  n_ddc_AVE, n_hc_AVE, n_hex_AVE, box_trans(cart,cart)
real :: n_ddc_AVE_SURF, n_hc_AVE_SURF, n_hex_AVE_SURF, n_ddc_AVE_BULK, n_hc_AVE_BULK, n_hex_AVE_BULK
real :: delta_AVE, delta_AVE_BULK, delta_AVE_SURF, esse_AVE, esse_AVE_BULK, esse_AVE_SURF, rog_AVE, rog_AVE_BULK, rog_AVE_SURF
real :: ze_AVE, ze_AVE_BULK, ze_AVE_SURF, hbdist, hbangle, tangle
real, allocatable :: pos(:,:), stat_nr_AVE(:), stat_nr_HB_AVE(:)
character*3 :: switch_cages, switch_hex, switch_r_cls, r_cls_W, switch_r_idx, switch_ffss, switch_hbck
character*4, allocatable :: sym(:)
character*100 :: rings_exe, natformat

! DFS stuff
integer :: ncr, mxvic, nat_cls, iat, jat, nnf, vol_count, voltot, ncrit, patch
integer :: crit=0 ! If we have something, we have an hexagon, so...
integer, allocatable :: dfs_color(:), volume_crit(:)

! Getting the - variable - box...
icell(1)=box_trans(1,1) ; icell(2)=box_trans(1,2) ; icell(3)=box_trans(1,3)
icell(4)=box_trans(2,1) ; icell(5)=box_trans(2,2) ; icell(6)=box_trans(2,3)
icell(7)=box_trans(3,1) ; icell(8)=box_trans(3,2) ; icell(9)=box_trans(3,3)

! Initialize the color array...
r_color(:)=0

! Write down an .xyz with the region we have selected
if (trim(adjustl(switch_r_idx)).eq.'no') then ! pick up those atoms within some z-slice
   open(unit=69, file='conf.xyz', status='unknown')
   open(unit=70, file='tmp.dat', status='unknown')
   nxyz=0
   kto(:)=0
   do i=1,r_ns
      do j=1,n_r_ws(i)
         if (pos(cart,list_r_ws(i,j)).ge.r_zmin.and.pos(cart,list_r_ws(i,j)).le.r_zmax) then
            nxyz=nxyz+1
            ! Index nxyz in conf.xyz corresponds to index list_r_ws(i,j) in the global .xtc
            ! Store this information for visualisation purposes
            kto(nxyz)=list_r_ws(i,j)
         endif
      enddo 
   enddo 
   write(70,*) nxyz
   write(70,'(4f20.10,1i10)') icell(1)*10.0, icell(5)*10.0, icell(9)*10.0, rcut*10.0, maxr_RINGS
   do i=1,r_ns
      do j=1,n_r_ws(i)
         if (pos(cart,list_r_ws(i,j)).ge.r_zmin.and.pos(cart,list_r_ws(i,j)).le.r_zmax) then
            !write(69,'(1a5,3f20.10)') sym(list_r_ws(i,j)), pos(:,list_r_ws(i,j))*10.0
            ! We write down everything as oxygen atoms!! 
            write(69,'(1a5,3f20.10)') "O", pos(:,list_r_ws(i,j))*10.0
         endif
      enddo
   enddo
   close(69)
   close(70)
else ! we have already read the indexes of the atoms we are interested in - typically some ice cluster...
     ! this number changes in time! TBF

   !! DEBUG
   !write(*,*) "HERE"
   !! END DEBUG

   open(unit=69, file='conf.xyz', status='unknown')
   open(unit=70, file='tmp.dat', status='unknown')
   nxyz=0
   kto(:)=0
   minz=1.0d30
   maxz=-1.0d30
   do i=1,C_size(counter+1)
      nxyz=nxyz+1
      ! Index nxyz in conf.xyz corresponds to index list_r_ws(i,j) in the global .xtc
      ! Store this information for visualisation purposes
      kto(nxyz)=C_idx(counter+1,i)+1

      !write(*,*) kto(nxyz)-1      

      ! color the atoms in the cluster...
      r_color(C_idx(counter+1,i)+1)=9
      ! Check whether this is a cluster at the surface or not
      if (pos(3,kto(nxyz)).lt.minz) minz=pos(3,kto(nxyz)) 
      if (pos(3,kto(nxyz)).gt.maxz) maxz=pos(3,kto(nxyz))
   enddo
   surfF=0
   ! Get the spatial extent of the clusters with respect to the z direction
   z_ext=(maxz-minz)
   ze_AVE=ze_AVE+z_ext

   !! JPCL 
   !open(unit=6251, file='jpcl_z_1', status='unknown', position="append")
   !open(unit=6252, file='jpcl_z_2', status='unknown', position="append")
   !open(unit=6253, file='jpcl_z_3', status='unknown', position="append")
   !open(unit=6254, file='jpcl_z_4', status='unknown', position="append")
   !!

   if (minz.le.thrS.or.maxz.ge.thrSS) then
      surfF=1 ! This is a surface cluster !
      nsurf=nsurf+1
      ze_AVE_SURF=ze_AVE_SURF+z_ext 
      !! JPCL stuff
      !if (mflag(counter+1).eq.1) then ! Surviving
      !   write(6251,*) counter, z_ext
      !endif !else
      !   write(6252,*) counter, z_ext ! All
      !!endif
      !!
   else
      nbulk=nbulk+1
      ze_AVE_BULK=ze_AVE_BULK+z_ext
      !! JPCL stuff
      !if (mflag(counter+1).eq.1) then ! this one made it through
      !   write(6253,*) counter, z_ext
      !endif !else
      !   write(6254,*) counter, z_ext
      !!endif
      !!
   endif 

   close(6251)
   close(6252)
   close(6253)
   close(6254)

   write(70,*) nxyz
   write(70,'(4f20.10,1i10)') icell(1)*10.0, icell(5)*10.0, icell(9)*10.0, rcut*10.0, maxr_RINGS
   do i=1,C_size(counter+1)
      write(69,'(1a5,3f20.10)') "O", pos(:,C_idx(counter+1,i)+1)*10.0
   enddo
   close(69)
   close(70)
endif

! Bash stuff...
! We just write each specie as oxygen atoms in any case...
!command="cat conf.xyz | sed 's/OW/O/'> temp && mv temp conf.xyz"
!call system(command)
command="cat conf.xyz >> tmp.dat ; mv tmp.dat conf.xyz"
call system(command)
command="n_xyz=`wc -l conf.xyz | awk '{print $1-2}'` ; cat rings.in_TEMPLATE | sed ""s/NAT/$n_xyz/"" > rings.in"
call system(command)
command="c1=`head -2 conf.xyz | tail -1 | awk '{print $1}'` ; cat rings.in | sed ""s/ICELL1/$c1/"" > tmp.dat"
call system(command)
command="c2=`head -2 conf.xyz | tail -1 | awk '{print $2}'` ; cat tmp.dat | sed ""s/ICELL2/$c2/"" > rings.in"
call system(command)
command="c3=`head -2 conf.xyz | tail -1 | awk '{print $3}'` ; cat rings.in | sed ""s/ICELL3/$c3/"" > tmp.dat"
call system(command)
command="rc=`head -2 conf.xyz | tail -1 | awk '{print $4}'` ; cat tmp.dat | sed ""s/RCUT/$rc/"" > rings.in"
call system(command)
command="mr=`head -2 conf.xyz | tail -1 | awk '{print $5}'` ; cat rings.in | sed ""s/MAXR/$mr/"" > tmp.dat"
call system(command)
command="mv tmp.dat rings.in"
call system(command)
command="mv conf.xyz data ; cp options_TEMPLATE options"
call system(command)
command="rm -r -f rings.out tmp rstat bonds Walltime rings.dat r3-5.dat r4-5.dat r5-5.dat r6-5.dat"
call system(command)
call system(rings_exe // "rings.in > log 2>&1")
!call system(rings_exe // "rings.in")

! DEBUG
!write(*,*) rings_exe
!command="touch diocane.log"
!call system(command)
!stop
!

! Allocate stuff for stat
allocate(stat_nr(maxr-2)) ! this guy contains the number of rings of each n-membered type
                          ! e.g. if maxr=9, n. of 3,4,5,6,7,8 or 9 membered rings
! HB
allocate(stat_nr_HB(maxr-2))

! Get all the files we need...
nmem=0
stat_nr(:)=0
stat_nr_HB(:)=0
do k=3,maxr
   nmem=nmem+1
   ! if non-primitive rings, substitute liste-5 with liste-1
   command="./rstat/liste-5/r"
   ! if non-primitive rings, substitute -5.dat with -1.dat
   rst="-5.dat"
   write(rst2,*) k
   rst2=trim(adjustl(rst2))
   fcommand=trim(command)//trim(rst2)//trim(rst)
   inquire(file=fcommand, exist=exist)
   if (exist) then
      command="cp ./rstat/liste-5/r"
      ! if non-primitive rings, substitute -5.dat with -1.dat
      rst="-5.dat ."
      write(rst2,*) k
      rst2=trim(adjustl(rst2))
      fcommand=trim(command)//trim(rst2)//trim(rst)
      call system(fcommand) 

      command="r"
      ! if non-primitive rings, substitute -5.dat with -1.dat
      rst="-5.dat"
      fcommand=trim(command)//trim(rst2)//trim(rst)
      open(unit=69, file=trim(adjustl(fcommand)), status='old')
      nl=0
      endf=0
      do
      read(69,*,iostat=endf)
        if (endf==-1) exit
        nl=nl+1
      enddo     
      rewind(69)
      stat_nr(nmem)=nl
      close(69)
   else
      write(99,*) "Achtung! At time ", time, "no ", k, "-membered rings were found!"
      stat_nr(nmem)=0
   endif  
enddo

! Allocate this rather cumbersome array...
allocate(stat_wr%stat_wr_size(maxr-2))
nmem=0
do k=3,maxr
   nmem=nmem+1
   allocate(stat_wr%stat_wr_size(nmem)%mrings(stat_nr(nmem),k)) ! Apparently no need to deallocate
   ! Get the atoms belonging to each ring size
   if (stat_nr(nmem).gt.0) then
      command="r"
      write(rst2,*) k
      rst2=trim(adjustl(rst2))
      ! if non-primitive rings, substitute -5.dat with -1.dat
      rst="-5.dat"
      fcommand=trim(command)//trim(rst2)//trim(rst)
      open(unit=69, file=trim(adjustl(fcommand)), status='old')
      do kr=1,stat_nr(nmem)
         ! Get the atoms involved in each ring... 
         read(69,*) stat_wr%stat_wr_size(nmem)%mrings(kr,:)
         ! Colors
         if (wcol.eq.k) then
            r_color(kto(stat_wr%stat_wr_size(nmem)%mrings(kr,:)))=k
            !write(*,*) kto(stat_wr%stat_wr_size(nmem)%mrings(kr,:))-1
         endif   
      enddo 
      close(69)
   endif
enddo

! We want to know the fraction of each n-membered category that are actually wholly hydrogen bonded
nmem=0
if (trim(adjustl(switch_hbck)).eq.'yes') then
   do k=3,maxr
      nmem=nmem+1
      if (stat_nr(nmem).gt.0) then
         do kr=1,stat_nr(nmem)
            !write(*,*) nmem, k, kr, kto(stat_wr%stat_wr_size(nmem)%mrings(kr,:))
            hbflag(:)=0 ! If this guy is HB, we should have hbflag=2 for each element!
            ! Loop over each oxygen
            do l=1,k
               ! Loop over every other oxygen
               do m=l+1,k
                  !! M-L
                  if (sym(kto(stat_wr%stat_wr_size(nmem)%mrings(kr,l))).eq."OR1".or.sym(kto(stat_wr%stat_wr_size(nmem)%mrings(kr,l))).eq."OR2".or.sym(kto(stat_wr%stat_wr_size(nmem)%mrings(kr,l))).eq."OR3".or.sym(kto(stat_wr%stat_wr_size(nmem)%mrings(kr,l))).eq."OR4") then ! l is metaldehide - or an oxygen beloging to an etheric group . It can only receive HB!
                     dummy=0
                  else ! water-water or water-OH TIP4P/Ice, OW HW1 HW2 MW
                     ! H1
                     xdf=pos(1,kto(stat_wr%stat_wr_size(nmem)%mrings(kr,l))+1)-pos(1,kto(stat_wr%stat_wr_size(nmem)%mrings(kr,m)))
                     ydf=pos(2,kto(stat_wr%stat_wr_size(nmem)%mrings(kr,l))+1)-pos(2,kto(stat_wr%stat_wr_size(nmem)%mrings(kr,m)))
                     zdf=pos(3,kto(stat_wr%stat_wr_size(nmem)%mrings(kr,l))+1)-pos(3,kto(stat_wr%stat_wr_size(nmem)%mrings(kr,m)))
                     call images(cart,0,1,1,icell,xdf,ydf,zdf)
                     rij(1)=xdf ; rij(2)=ydf ; rij(3)=zdf; rij_M=sqrt(rij(1)**2.0+rij(2)**2.0+rij(3)**2.0)
!                     ! DEBUG
!                     if (kto(stat_wr%stat_wr_size(nmem)%mrings(kr,l))-1.eq.5008.and.kr.eq.13.and.k.eq.6) then
!                        write(*,*) k, kr, kto(stat_wr%stat_wr_size(nmem)%mrings(kr,l)), kto(stat_wr%stat_wr_size(nmem)%mrings(kr,m))-1, rij_M
!                        !write(*,*) k, kr, kto(stat_wr%stat_wr_size(nmem)%mrings(kr,:))-1
!                     endif
                     ! END DEBUG
                     if (rij_M.lt.hbdist) then
                     ! Check the O-H-O angle
                     v1(:)=pos(:,kto(stat_wr%stat_wr_size(nmem)%mrings(kr,l))+1)-pos(:,kto(stat_wr%stat_wr_size(nmem)%mrings(kr,l)))
                     v1m=sqrt(v1(1)**2.0+v1(2)**2.0+v1(3)**2.0)
                     v2(:)=pos(:,kto(stat_wr%stat_wr_size(nmem)%mrings(kr,m)))-pos(:,kto(stat_wr%stat_wr_size(nmem)%mrings(kr,l))+1)
                     v2m=sqrt(v2(1)**2.0+v2(2)**2.0+v2(3)**2.0)
                     tangle=acos(((v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3))/(v1m*v2m)))*rad2deg
                        if (abs(tangle).lt.hbangle) then
                           hbflag(l)=hbflag(l)+1
                           hbflag(m)=hbflag(m)+1
                        endif
                     endif
                     ! H2
                     if (sym(kto(stat_wr%stat_wr_size(nmem)%mrings(kr,l))).ne."O3") then ! if .eq. it means that's a OH...
                        xdf=pos(1,kto(stat_wr%stat_wr_size(nmem)%mrings(kr,l))+2)-pos(1,kto(stat_wr%stat_wr_size(nmem)%mrings(kr,m)))
                        ydf=pos(2,kto(stat_wr%stat_wr_size(nmem)%mrings(kr,l))+2)-pos(2,kto(stat_wr%stat_wr_size(nmem)%mrings(kr,m)))
                        zdf=pos(3,kto(stat_wr%stat_wr_size(nmem)%mrings(kr,l))+2)-pos(3,kto(stat_wr%stat_wr_size(nmem)%mrings(kr,m)))
                        call images(cart,0,1,1,icell,xdf,ydf,zdf)
                        rij(1)=xdf ; rij(2)=ydf ; rij(3)=zdf; rij_M=sqrt(rij(1)**2.0+rij(2)**2.0+rij(3)**2.0)
                        if (rij_M.lt.hbdist) then
                        ! Check the O-H-O angle
                        v1(:)=pos(:,kto(stat_wr%stat_wr_size(nmem)%mrings(kr,l))+2)-pos(:,kto(stat_wr%stat_wr_size(nmem)%mrings(kr,l)))
                        v1m=sqrt(v1(1)**2.0+v1(2)**2.0+v1(3)**2.0)
                        v2(:)=pos(:,kto(stat_wr%stat_wr_size(nmem)%mrings(kr,m)))-pos(:,kto(stat_wr%stat_wr_size(nmem)%mrings(kr,l))+2)
                        v2m=sqrt(v2(1)**2.0+v2(2)**2.0+v2(3)**2.0)
                        tangle=acos(((v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3))/(v1m*v2m)))*rad2deg
                           if (abs(tangle).lt.hbangle) then
                              hbflag(l)=hbflag(l)+1
                              hbflag(m)=hbflag(m)+1
                           endif
                        endif
                     endif
                  endif
                  !! L-M
                  if (sym(kto(stat_wr%stat_wr_size(nmem)%mrings(kr,m))).eq."OR1".or.sym(kto(stat_wr%stat_wr_size(nmem)%mrings(kr,m))).eq."OR2".or.sym(kto(stat_wr%stat_wr_size(nmem)%mrings(kr,m))).eq."OR3".or.sym(kto(stat_wr%stat_wr_size(nmem)%mrings(kr,m))).eq."OR4") then ! l is metaldehide - or an oxygen beloging to an etheric group
                     dummy=dummy
                  else ! water-water or water-OH
                     ! H1
                     xdf=pos(1,kto(stat_wr%stat_wr_size(nmem)%mrings(kr,m))+1)-pos(1,kto(stat_wr%stat_wr_size(nmem)%mrings(kr,l)))
                     ydf=pos(2,kto(stat_wr%stat_wr_size(nmem)%mrings(kr,m))+1)-pos(2,kto(stat_wr%stat_wr_size(nmem)%mrings(kr,l)))
                     zdf=pos(3,kto(stat_wr%stat_wr_size(nmem)%mrings(kr,m))+1)-pos(3,kto(stat_wr%stat_wr_size(nmem)%mrings(kr,l)))
                     call images(cart,0,1,1,icell,xdf,ydf,zdf)
                     rij(1)=xdf ; rij(2)=ydf ; rij(3)=zdf; rij_M=sqrt(rij(1)**2.0+rij(2)**2.0+rij(3)**2.0)
                     if (rij_M.lt.hbdist) then
                     ! Check the O-H-O angle
                     v1(:)=pos(:,kto(stat_wr%stat_wr_size(nmem)%mrings(kr,m))+1)-pos(:,kto(stat_wr%stat_wr_size(nmem)%mrings(kr,m)))
                     v1m=sqrt(v1(1)**2.0+v1(2)**2.0+v1(3)**2.0)
                     v2(:)=pos(:,kto(stat_wr%stat_wr_size(nmem)%mrings(kr,l)))-pos(:,kto(stat_wr%stat_wr_size(nmem)%mrings(kr,m))+1)
                     v2m=sqrt(v2(1)**2.0+v2(2)**2.0+v2(3)**2.0)
                     tangle=acos(((v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3))/(v1m*v2m)))*rad2deg
                        if (abs(tangle).lt.hbangle) then
                           hbflag(l)=hbflag(l)+1
                           hbflag(m)=hbflag(m)+1
                        endif
                     endif
                     ! H2
                     if (sym(kto(stat_wr%stat_wr_size(nmem)%mrings(kr,m))).ne."O3") then ! if .eq. it means that's a OH...
                        xdf=pos(1,kto(stat_wr%stat_wr_size(nmem)%mrings(kr,m))+2)-pos(1,kto(stat_wr%stat_wr_size(nmem)%mrings(kr,l)))
                        ydf=pos(2,kto(stat_wr%stat_wr_size(nmem)%mrings(kr,m))+2)-pos(2,kto(stat_wr%stat_wr_size(nmem)%mrings(kr,l)))
                        zdf=pos(3,kto(stat_wr%stat_wr_size(nmem)%mrings(kr,m))+2)-pos(3,kto(stat_wr%stat_wr_size(nmem)%mrings(kr,l)))
                        call images(cart,0,1,1,icell,xdf,ydf,zdf)
                        rij(1)=xdf ; rij(2)=ydf ; rij(3)=zdf; rij_M=sqrt(rij(1)**2.0+rij(2)**2.0+rij(3)**2.0)
                        if (rij_M.lt.hbdist) then
                        ! Check the O-H-O angle
                        v1(:)=pos(:,kto(stat_wr%stat_wr_size(nmem)%mrings(kr,m))+2)-pos(:,kto(stat_wr%stat_wr_size(nmem)%mrings(kr,m)))
                        v1m=sqrt(v1(1)**2.0+v1(2)**2.0+v1(3)**2.0)
                        v2(:)=pos(:,kto(stat_wr%stat_wr_size(nmem)%mrings(kr,l)))-pos(:,kto(stat_wr%stat_wr_size(nmem)%mrings(kr,m))+2)
                        v2m=sqrt(v2(1)**2.0+v2(2)**2.0+v2(3)**2.0)
                        tangle=acos(((v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3))/(v1m*v2m)))*rad2deg
                           if (abs(tangle).lt.hbangle) then
                              hbflag(l)=hbflag(l)+1
                              hbflag(m)=hbflag(m)+1
                           endif
                        endif
                     endif
                  endif     
               enddo
            enddo 
            !write(*,*) kto(stat_wr%stat_wr_size(nmem)%mrings(kr,:))-1
            !write(*,*) (hbflag(l), l=1,k)
            ! HERE !
            !! stat_nr contains the number of rings / n-member. if all hbflags = 2 , this is HB - DONE
            !! so we count a stat_nr_HB, average that as well, and output that as well - DONE
            !! also, we should fix the color - DONE 
            !! also, we should fix OW and O3, different options! - DONE. the O of MDHE will have to be implemented separately...
            !! 
            if (sum(hbflag).eq.(2*k)) then
               !write(*,*) k, sum(hbflag)
               stat_nr_HB(nmem)=stat_nr_HB(nmem)+1
               ! fix the color!
               if (wcol.eq.k) then
                  r_color(kto(stat_wr%stat_wr_size(nmem)%mrings(kr,:)))=k*100 
               endif 
               !if (k.eq.6) then
               !   write(*,*) kto(stat_wr%stat_wr_size(nmem)%mrings(kr,:))-1
               !endif
            endif
         enddo
      endif      
   enddo
endif

! TBF 
!stop
! TBF

! Hexagonal rings and cages...

!!if (trim(adjustl(switch_cages)).eq.'yes'.and.trim(adjustl(switch_hex)).eq.'no') then
!!   write(99,*) "You need to compute hexagonal rings to get the cages... i.e. CAGES=yes and HEX=yes"
!!   ! Is this actually true?! NO
!!   stop
!!endif

if (trim(adjustl(switch_cages)).eq.'yes'.or.trim(adjustl(switch_hex)).eq.'yes') then
    ! Do some checks...
    if (maxr.lt.6) then
       write(99,*) "You need to get 6-membered rings in order to look at hexagons and/or DDCs and HCs..."
       stop
    endif 
   
   ! Six-rings stuff, for hexagons/cages. Keep it separate from "normal" rings statistics...
   ! if non-primitive rings, substitute -5.dat with -1.dat
   if (stat_nr(4)==0) then ! No hexagonal rings!?
      n_ddc=0
      n_hc=0
      goto 656
   endif
   open(unit=69, file='r6-5.dat', status='old')
   nl=0
   endf=0
   do 
   read(69,*,iostat=endf)
     if (endf==-1) exit
     nl=nl+1
   enddo      
   ! N. of rings = nl
   allocate(w_rings(nl,6),r_array(nl*6),p_rings(6,nl,6),r_nper(6)) ! Hexagonal rings only, in here...
   rewind(69)
   do k=1,nl
      ! Get the atoms involved in each ring... 
      read(69,*) w_rings(k,:)
      do l=1,6
         id=((k-1)*6)+l
         r_array(id)=w_rings(k,l)
      enddo
   enddo

   if (trim(adjustl(switch_hex)).eq.'yes') then
      n_hex=0
      !! DEBUG
      !write(*,*) counter, nl
      !! END DEBUG 
      do k=1,nl
         ! Loop over the 6 possible triplets:
         r_flag=0
         do tri=1,6
            ti=tri
            if ((tri+1).eq.6) then
               tj=6
            else
               tj=mod((tri+1),6)
            endif
            if ((tri+2).eq.6) then
               tk=6
            else
               tk=mod((tri+2),6)
            endif
            ! Get the angles for each triplet
            posi(:)=pos(:,kto(w_rings(k,ti)))
            posj(:)=pos(:,kto(w_rings(k,tj)))
            posk(:)=pos(:,kto(w_rings(k,tk)))
            ! For some reason, some rings just don't make ANY sense. Discard them on the basis of angles AND distances...
            xdf=posi(1)-posj(1) ; ydf=posi(2)-posj(2) ; zdf=posi(3)-posj(3)
            call images(cart,0,1,1,icell,xdf,ydf,zdf)
            rij(1)=xdf ; rij(2)=ydf ; rij(3)=zdf; rij_M=sqrt(rij(1)**2.0+rij(2)**2.0+rij(3)**2.0)
            xdf=posk(1)-posj(1) ; ydf=posk(2)-posj(2) ; zdf=posk(3)-posj(3)
            call images(cart,0,1,1,icell,xdf,ydf,zdf)
            rkj(1)=xdf ; rkj(2)=ydf ; rkj(3)=zdf; rkj_M=sqrt(rkj(1)**2.0+rkj(2)**2.0+rkj(3)**2.0)
            hex_angle=acos(((rij(1)*rkj(1)+rij(2)*rkj(2)+rij(3)*rkj(3))/(rij_M*rkj_M)))*rad2deg
            if ((ref_angle-hex_angle).gt.a_thr) then
               r_flag=1 !; exit # this one saves ome time...
            endif
         enddo
         if (r_flag.eq.0) then ! We do have a proper hexagonal ring
            n_hex=n_hex+1 

            !! DEBUG
            !write(*,*) counter, nl, n_hex
            !! END DEBUG 
 
            if (wcol.eq.0.or.wcol.eq.1) then ! TBF
               r_color(kto(w_rings(k,:)))=10
            endif
         endif
      enddo
   endif

   n_hex_AVE=n_hex_AVE+real(n_hex)
   ! Time [ps] | N. of 6-membered rings | N. of proper hexagonal rings
   if (trim(adjustl(switch_hex)).eq.'yes') then
      open(unit=108, file='hin_structure.out.rings.hex', status='unknown', position='append')
      write(108,'(1E10.4,2i15)') time, nl, n_hex
      close(108)
   endif

   ! Cluster hexagonal rings, e.g. to find the largest patch of hexagonal rings sitting on top of the surface
   if (trim(adjustl(switch_r_cls)).eq.'yes') then
      if (trim(adjustl(switch_hex)).ne.'yes') then
         write(99,*) "You can only cluster regular hexagonal rings at the moment..."
         stop
      endif
      if (trim(adjustl(r_cls_W)).eq.'CLA') then
         call clath_cages(stat_wr,stat_nr,time,nat,natformat,kto)
         
      else if (trim(adjustl(r_cls_W)).ne.'SIX') then
         write(99,*) "Sorry mate, I can do only six membered rings at the moment..."
         stop
      else
      !
      allocate(cr_list(nat),lwho(nat,nat),dfs_color(nat),volume_crit(nat))
      ncr=0
      cr_list(:)=0
      do i=1,nat
         if (r_color(i).eq.10) then
            ncr=ncr+1
            cr_list(ncr)=i
         endif
      enddo
      ! allocate the adjecency list and related stuff..
      mxvic=6 ! at most, atoms involved in an hexagonal ring can have 6 atoms
      nat_cls=nl*6 ! the biggest cluster possible involves all the atoms involved with hex rings, nl*6
      allocate(graph_solid_connect(ncr,mxvic),followgraph(ncr),volume(nat_cls))
      allocate(neigh(ncr),predecessor(ncr))
      ! fill the adjacency list...
      neigh(:)=0
      lwho(:,:)=0
      graph_solid_connect(:,:)=0
      predecessor(:)=0

      do j=1,ncr
         do k=1,ncr
            if (j.ne.k) then
               iat=cr_list(j)
               jat=cr_list(k)
               posj(:)=pos(:,jat) ; posi(:)=pos(:,iat)
               call nn (posj,posi,icell,rcut,cknn,xdf,ydf,zdf,dist)
               if (cknn) then 
                  neigh(j)=neigh(j)+1
                  graph_solid_connect(j,neigh(j))=k
               endif
            endif
         enddo
      enddo

      count_cls=0
      volume(:)=0
      vol_count=0
      followgraph(:)=0

      do j=1,ncr
         black=0
         if(followgraph(j).eq.0)then
            count_cls=count_cls+1
            followgraph(j)=explore(j)
            volume(count_cls)=black
         end if
      end do

      nnf=count_cls

      ! number of clusters (not yet filtered by crit) = nnf  
      ! volume of each cluster = volume(j) up to count
      ! list of atoms indexes that constitute the clusters = lwho

      count_cls=0
      dfs_color(:)=0
      voltot=0
      volume_crit(:)=0

      do j=1,nnf
         if (volume(j).gt.crit) then
            do k=1,volume(j)
               dfs_color(lwho(j,k))=j
            enddo
            count_cls=count_cls+1
            voltot=voltot+volume(j)
            volume_crit(count_cls)=volume(j)
         endif
      enddo

     ncrit=count_cls

     ! number of clusters greater than crit = ncrit

     ! sort the volumes: the biggest is the surface itself
     call sort2(volume_crit,count_cls)

     ! Here is the number of atoms within the largest hexagonal patch
     patch=volume_crit(1)
     ! Write down patch statistics: time, n. of atoms in the biggest patch, n. of atoms involved in hex rings as a whole
     write(207,"(1E10.4,2i10)") time, patch, ncr

     ! Write down the colors for VMD
     write(208,"("//adjustl(natformat)//"i10)") (dfs_color(k), k=1,nat)

!     ! output by columns:
!     ! col0 = time [ps]
!     ! col1 = number of crystalline particles = ncr
!     ! col2 = number of crystalline particles involved in crystalline clusters of volume > crit = voltot
!     ! col3 = number of crystalline clusters of volume > crit = ncrit
!     ! col4 = volume of the biggest cluster (the surface itself) = volume_crit(1)
!     ! col5-> up to ncrit = volumes of each crit cluster
!     write(300,'(1f12.6,3i10,<count>i10)') dble(i)*tframe, ncr, voltot, ncrit, (volume_crit(j), j=1,count)
!
!     ! color
!     write(301,'(<nat>i10)') (color(j), j=1,nat)
!     !write(301,'(<nat>i10)') (oflag(j), j=1,nat)

!!



      deallocate(graph_solid_connect,followgraph,volume,neigh,predecessor,lwho,cr_list,dfs_color,volume_crit)
      endif
      !
   endif
   ! 

   if (trim(adjustl(switch_cages)).eq.'yes') then
         
      ! Find DDCs cages
      ! 1. For every 1<k<6, a minimum of three other hexagonal rings pass through mk.
      ! 2. For every triplet Tk = (mk, mk@1, mk@2), there is at least one hexagonal ring
      ! other than R0 that passes through mk, mk@1 and mk@2. Here a@b = (a+b).mod.6
      ! There are 6 triplets possible for a six-membered ring
      ! 3. The three top peripheral rings, R1,R3 and R5, and the three bottom
      ! peripheral rings, R2,R4 and R6, each have one water molecule in
      ! common. 
      n_ddc=0
      do k=1,nl 
         r_flag=0
         do l=1,6
            nsix=count(r_array.eq.w_rings(k,l))
            if (nsix.lt.4) r_flag=1
         enddo 
         if (r_flag.eq.0) then
            ! 1. is fulfilled. Now, get the triplets for 2.
            r_flag3=0
            do l=1,6 
               tr(1)=w_rings(k,l) 
               if (l+1.ne.6) then
                  tr(2)=w_rings(k,mod((l+1),6))
               else
                  tr(2)=w_rings(k,6)
               endif
               if (l+2.ne.6) then
                  tr(3)=w_rings(k,mod((l+2),6))
               else
                  tr(3)=w_rings(k,6)
               endif
               ! For every triplet, check whether there is at least one OTHER
               ! hexagonal ring involving the triplet
               r_flag2=1
               nper=0
               do m=1,nl
                  if (m.eq.k) cycle
                  tr6(:)=w_rings(m,:)
                  if ((count(tr6.eq.tr(1))+count(tr6.eq.tr(2))+count(tr6.eq.tr(3))).eq.3) then
                     r_flag2=0
                     ! Store the information about the peripheral rings
                     ! p_rings(6,nl,6) rj,ri,rk. ri-> more than one ring
                     ! passing through the triplet (nper).  rj -> which triplet (l)
                     ! . rk -> members of the peripheral ring
                     nper=nper+1
                     r_nper(l)=nper ! How many peripheral rings for triplet l
                     p_rings(l,nper,:)=tr6(:) ! Keep track of the atoms involved in each p-ring...
                  endif
               enddo
               if (r_flag2.ne.0) r_flag3=1
            enddo ! Loop over triplets
            if (r_flag3.eq.0) then
               ! 2. is fulfilled, now check 3.
               ! Odd peripherals
               do per1=1,r_nper(1) ! loop over the p-rings of T1
                  do per3=1,r_nper(3) ! loop over the p-rings of T3
                     do kper135=1,6
                        r13=count(p_rings(3,per3,:).eq.p_rings(1,per1,kper135))
                        if (r13.gt.0) then ! R1 and R3 have one molecule (kper) in common
                           do per5=1,r_nper(5)
                              r15=count(p_rings(5,per5,:).eq.p_rings(1,per1,kper135))
                              if (r15.gt.0) then ! R1, R3 and R5 have one molecule in common
                                 ! Even peripherals
                                 do per2=1,r_nper(2) ! loop over the p-rings of T2
                                    do per4=1,r_nper(4) ! loop over the p-rings of T4
                                       do kper246=1,6
                                          r24=count(p_rings(4,per4,:).eq.p_rings(2,per2,kper246))
                                          if (r24.gt.0) then ! R2 and R4 have one molecule (kper) in common
                                             do per6=1,r_nper(6)
                                                r26=count(p_rings(6,per6,:).eq.p_rings(2,per2,kper246))
                                                if (r26.gt.0) then ! R2, R4 and R6 have one molecule in common
                                                   ! Now, check that r1|r3, r1|r5 , r3|r5 -> ! three molecules in ! common
                                                   ! As well as r2|r4, r2|r6 ! , r4|r6     
                                                   ckr=0
                                                   do ck=1,6  
                                                      ckr=ckr+count(p_rings(3,per3,:).eq.p_rings(1,per1,ck))
                                                   enddo 
                                                   if (ckr.ne.3) r_flag3=1      
                                                   ckr=0
                                                   do ck=1,6
                                                      ckr=ckr+count(p_rings(5,per5,:).eq.p_rings(1,per1,ck))
                                                   enddo
                                                   if (ckr.ne.3) r_flag3=1
                                                   ckr=0
                                                   do ck=1,6
                                                      ckr=ckr+count(p_rings(5,per5,:).eq.p_rings(3,per3,ck))
                                                   enddo
                                                   if (ckr.ne.3) r_flag3=1
                                                   ckr=0
                                                   do ck=1,6
                                                      ckr=ckr+count(p_rings(4,per4,:).eq.p_rings(2,per2,ck))
                                                   enddo
                                                   if (ckr.ne.3) r_flag3=1
                                                   ckr=0
                                                   do ck=1,6
                                                      ckr=ckr+count(p_rings(6,per6,:).eq.p_rings(2,per2,ck))
                                                   enddo
                                                   if (ckr.ne.3) r_flag3=1
                                                   ckr=0
                                                   do ck=1,6
                                                      ckr=ckr+count(p_rings(6,per6,:).eq.p_rings(4,per4,ck))
                                                   enddo
                                                   if (ckr.ne.3) r_flag3=1
                                                   if (r_flag3.eq.0) then ! we finally have a proper DDC cage
                                                       n_ddc=n_ddc+1 
                                                       if (wcol.eq.1) then
                                                               r_color(kto(w_rings(k,:)))=15 
                                                          r_color(kto(p_rings(1,per1,:)))=15 
                                                          r_color(kto(p_rings(3,per3,:)))=15
                                                          r_color(kto(p_rings(5,per5,:)))=15
                                                          r_color(kto(p_rings(2,per2,:)))=15
                                                          r_color(kto(p_rings(4,per4,:)))=15
                                                          r_color(kto(p_rings(6,per6,:)))=15
                                                       endif
                                                   endif
                                                endif
                                             enddo
                                          endif
                                       enddo
                                    enddo
                                 enddo
                              ! end Even
                              endif
                           enddo
                        endif
                     enddo   
                  enddo
               enddo
            endif
         endif
      enddo
                  
      ! HCs rings...
      ! Find HCs cages
      ! 1. R1 | R2 == 0
      ! 2. There exists 1.le.k.le.6 so that mk is a neighbor of l1 or l2
      ! as defined based on the distance criterion.
      ! 3. If mk is a neighbour of l1, mk2 and mk4 must be neighbours of l3 and l5 (or l5 and l3), respectively.
      ! Adjusting the algorithm to the case of mk being a neighbor of l2 is straightforward.
      n_hc=0

!! DEBUG
!write(*,*) "DEBUG"
!do i=1,nl
!   write(*,*) i, nl, kto(w_rings(i,:))-1
!enddo
!! END DEBUG


      do l=1,nl
         do m=1,nl
            if (l.eq.m) cycle
            r_flag=0
            do k=1,6
               if (count(w_rings(l,:).eq.w_rings(m,k)).gt.0) r_flag=1
            enddo
            if (r_flag.eq.0) then ! 1. is fulfilled
               r_flag2=1
               do k=1,6 
                  ! mk nn l1
                  ! DEBUG
!                  if (kto(w_rings(l,1)).eq.0) then
                  !if (kto(w_rings(m,k)).eq.0) then
                  !    !write(*,*) "ACHTUNG! Some silly stuff is happening for this particular configuration..."
                  !    write(99,*) "ACHTUNG! Some silly stuff is happening for this particular configuration..."
                  !    n_ddc=-1
                  !    n_hc=-1
                  !    goto 656 ! shit happened!
                  ! endif 
                   !  write(*,*) "DEBUG"
                   !  do i=1,nl
                   !     write(*,*) i, nl, w_rings(i,:)
                   !  enddo
                     !write(*,*) nl, w_rings(m,k), kto(w_rings(m,k))
                     !stop 
                  !endif 
                  ! END DEBUG
                  posj(:)=pos(:,kto(w_rings(l,1))) ; posi(:)=pos(:,kto(w_rings(m,k)))
                  call nn (posj,posi,icell,rcut,cknn,xdf,ydf,zdf,dist)
                  if (cknn) then ! mk is a neighbor of l1 | 2. is fulfilled
                     ! Check whether mk+2 nn l3 and mk+4 nn l5 OR mk+2 nn l5 and mk+4 nn l3
                     ! mk+2 nn l3
                     posj(:)=pos(:,kto(w_rings(l,3)))
                     if (k+2.ne.6) then
                        posi(:)=pos(:,kto(w_rings(m,mod((k+2),6)))) 
                     else
                        posi(:)=pos(:,kto(w_rings(m,6)))  
                     endif
                     call nn (posj,posi,icell,rcut,cknn,xdf,ydf,zdf,dist)
                     if (cknn) then 
                        ! mk+4 nn l5
                        posj(:)=pos(:,kto(w_rings(l,5)))
                        if (k+4.ne.6) then
                           posi(:)=pos(:,kto(w_rings(m,mod((k+4),6))))
                        else
                           posi(:)=pos(:,kto(w_rings(m,6)))
                        endif
                        call nn (posj,posi,icell,rcut,cknn,xdf,ydf,zdf,dist)
                        if (cknn) then
                           r_flag2=0 ! 3. is fulfilled 
                        endif 
                     endif
                     ! OR...
                     ! mk+2 nn l5
                     posj(:)=pos(:,kto(w_rings(l,5)))
                     if (k+2.ne.6) then
                        posi(:)=pos(:,kto(w_rings(m,mod((k+2),6))))
                     else
                        posi(:)=pos(:,kto(w_rings(m,6)))
                     endif
                     call nn (posj,posi,icell,rcut,cknn,xdf,ydf,zdf,dist)
                     if (cknn) then
                        ! mk+4 nn l3
                        posj(:)=pos(:,kto(w_rings(l,3)))
                        if (k+4.ne.6) then
                           posi(:)=pos(:,kto(w_rings(m,mod((k+4),6))))
                        else
                           posi(:)=pos(:,kto(w_rings(m,6)))
                        endif
                        call nn (posj,posi,icell,rcut,cknn,xdf,ydf,zdf,dist)
                        if (cknn) then
                           r_flag2=0 ! 3. is fulfilled 
                        endif
                     endif
                  endif
                  ! mk nn l2
                  posj(:)=pos(:,kto(w_rings(l,2))) ; posi(:)=pos(:,kto(w_rings(m,k)))
                  call nn (posj,posi,icell,rcut,cknn,xdf,ydf,zdf,dist)
                  if (cknn) then ! mk is a neighbor of l2 | 2. is fulfilled
                     ! Check whether mk+2 nn l4 and mk+4 nn l6 OR mk+2 nn l6 and mk+4 nn l2
                     ! mk+2 nn l4
                     posj(:)=pos(:,kto(w_rings(l,4)))
                     if (k+2.ne.6) then
                        posi(:)=pos(:,kto(w_rings(m,mod((k+2),6))))
                     else
                        posi(:)=pos(:,kto(w_rings(m,6)))
                     endif
                     call nn (posj,posi,icell,rcut,cknn,xdf,ydf,zdf,dist)
                     if (cknn) then
                        ! mk+4 nn l6
                        posj(:)=pos(:,kto(w_rings(l,6)))
                        if (k+4.ne.6) then
                           posi(:)=pos(:,kto(w_rings(m,mod((k+4),6))))
                        else
                           posi(:)=pos(:,kto(w_rings(m,6)))
                        endif
                        call nn (posj,posi,icell,rcut,cknn,xdf,ydf,zdf,dist)
                        if (cknn) then
                           r_flag2=0 ! 3. is fulfilled 
                        endif
                     endif
                     ! OR...
                     ! mk+2 nn l6
                     posj(:)=pos(:,kto(w_rings(l,6)))
                     if (k+2.ne.6) then
                        posi(:)=pos(:,kto(w_rings(m,mod((k+2),6))))
                     else
                        posi(:)=pos(:,kto(w_rings(m,6)))
                     endif
                     call nn (posj,posi,icell,rcut,cknn,xdf,ydf,zdf,dist)
                     if (cknn) then
                        ! mk+4 nn l4
                        posj(:)=pos(:,kto(w_rings(l,4)))
                        if (k+4.ne.6) then
                           posi(:)=pos(:,kto(w_rings(m,mod((k+4),6))))
                        else
                           posi(:)=pos(:,kto(w_rings(m,6)))
                        endif
                        call nn (posj,posi,icell,rcut,cknn,xdf,ydf,zdf,dist)
                        if (cknn) then
                           r_flag2=0 ! 3. is fulfilled 
                        endif
                     endif
                  endif   
               enddo
               if (r_flag2.eq.0) then ! ! we finally have a proper HC cage
                  n_hc=n_hc+1
                  if (wcol.eq.1) then
                     r_color(kto(w_rings(l,:)))=20
                     r_color(kto(w_rings(m,:)))=20
                  endif
               endif
            endif
         enddo
      enddo
      deallocate(w_rings,r_array,p_rings,r_nper)
      close(69)

      !! JPCL 
      !open(unit=4251, file='jpcl_cages_1', status='unknown', position="append")
      !open(unit=4252, file='jpcl_cages_2', status='unknown', position="append")
      !open(unit=4253, file='jpcl_cages_3', status='unknown', position="append")
      !open(unit=4254, file='jpcl_cages_4', status='unknown', position="append")
      !open(unit=4255, file='jpcl_cages_5', status='unknown', position="append")
      !open(unit=4256, file='jpcl_cages_6', status='unknown', position="append")
      !open(unit=4257, file='jpcl_cages_7', status='unknown', position="append")
      !open(unit=4258, file='jpcl_cages_8', status='unknown', position="append")
      !!

      
      n_ddc_AVE=n_ddc_AVE+real(n_ddc)
      n_hc_AVE=n_hc_AVE+real(n_hc)
      if (surfF.eq.0) then ! within the bulk
         n_ddc_AVE_BULK=n_ddc_AVE_BULK+real(n_ddc)
         n_hc_AVE_BULK=n_hc_AVE_BULK+real(n_hc)
         !! JPCL stuff
         !if (mflag(counter+1).eq.1) then ! this one made it through
         !   write(4251,*) counter, n_ddc
         !   write(4252,*) counter, n_hc
         !endif !else
         !   write(4253,*) counter, n_ddc
         !   write(4254,*) counter, n_hc
         !!endif
         !!
      else ! at the surface
         n_ddc_AVE_SURF=n_ddc_AVE_SURF+real(n_ddc)
         n_hc_AVE_SURF=n_hc_AVE_SURF+real(n_hc)
         !! JPCL stuff
         !if (mflag(counter+1).eq.1) then ! this one made it through
         !   write(4255,*) counter, n_ddc
         !   write(4256,*) counter, n_hc
         !endif !else
         !   write(4257,*) counter, n_ddc
         !   write(4258,*) counter, n_hc
         !!endif
         !!
      endif 

      !! JPCL
      !close(4251)
      !close(4252)
      !close(4253)
      !close(4254)
      !close(4255)
      !close(4256)
      !close(4257)
      !close(4258)
      !!


      ! Time [ps] | N. of 6-membered rings | N. of DDC cages | N. of HC cages     
656   write(103,'(1E10.4,3i15)') time, nl, n_ddc, n_hc
   endif ! cages
endif ! hex/cages

if (adjustl(trim(switch_ffss)).eq.'yes') then
   write(*,*) "Currently not implemented - missing libs on MacOS..."
   write(99,*) "Currently not implemented - missing libs on MacOS..."
   !!!! Get the asphericity (\Delta) and the shape (S) parameters for the nucleus. 
   !!!
   !!!! com
   !!!
   ! Center of mass, all masses=1
   rcm(:)=0.0
   pol1(:)=0.0
   pol2(:)=0.0
   cell_ORTO(1)=icell(1) ; cell_ORTO(2)=icell(5) ; cell_ORTO(3)=icell(9)  
   
   !! DEBUG
   !write(*,*) "HERE", counter, C_size(counter+1)
   !! END DEBUG
   
   do k=1,C_size(counter+1) ! loop over all the vertexes k of the void j
      do l=1,cart
         ! back to com...
         pol1(l)=pol1(l)+(cell_ORTO(l)*cos((pos(l,kto(k))/cell_ORTO(l))*2.0*pi))
         pol2(l)=pol2(l)+(cell_ORTO(l)*sin((pos(l,kto(k))/cell_ORTO(l))*2.0*pi))
      enddo
   enddo
   do l=1,cart
      pol1(l)=pol1(l)/real(C_size(counter+1))
      pol2(l)=pol2(l)/real(C_size(counter+1))
      rcm(l)=((atan2(-pol2(l),-pol1(l)))+pi)/(2.0*pi)*cell_ORTO(l)
   enddo
   
   !! DEBUG
   !write(6996,*) counter, rcm(cart)*10.0 ! z component in angs...
   
   tin(:,:)=0.0
   do k=1,C_size(counter+1)
      xdf=pos(1,kto(k))-rcm(1)
      ydf=pos(2,kto(k))-rcm(2)
      zdf=pos(3,kto(k))-rcm(3)
      call images(cart,0,1,1,icell,xdf,ydf,zdf)
      tin(1,1)=tin(1,1)+(xdf*xdf)
      tin(1,2)=tin(1,2)+(xdf*ydf)
      tin(1,3)=tin(1,3)+(xdf*zdf)
      tin(2,1)=tin(2,1)+(ydf*xdf)
      tin(2,2)=tin(2,2)+(ydf*ydf)
      tin(2,3)=tin(2,3)+(ydf*zdf)
      tin(3,1)=tin(3,1)+(zdf*xdf)
      tin(3,2)=tin(3,2)+(zdf*ydf)
      tin(3,3)=tin(3,3)+(zdf*zdf)
   enddo
   tin(:,:)=tin(:,:)/real(2.0*((real(C_size(counter+1)))**2.0))
   mtemp(:,:)=dble(tin(:,:))
   !!!!! DEBUG
   !!write(*,*) counter
   !!write(*,*) tin(1,:)
   !!write(*,*) tin(2,:)
   !!write(*,*) tin(3,:)
   !!write(*,*) mtemp(1,:)
   !!write(*,*) mtemp(2,:)
   !!write(*,*) mtemp(3,:)
   !!!stop
   !!!!! END
   !
   lwork=1000
   allocate(work(lwork))
   ! Commenting just this guy - missing libs on MacOS... TBF
   !call dsyev('N','U',cart,mtemp,cart,eigen,work,lwork,info)
   deallocate(work)
   ! DEBUG
   !write(*,*) "EIGEN", eigen(:) 
   ! END DEBUG
   trt=eigen(1)+eigen(2)+eigen(3)
   rog=dsqrt(trt)
   trt2=trt**2.0
   lambda=trt/3.0
   delta=(3.0/(2.0*trt2))*((eigen(1)-lambda)**2.0d0+(eigen(2)-lambda)**2.0d0+(eigen(3)-lambda)**2.0d0)
   esse=(27.0/(trt**3.0))*((eigen(1)-lambda)*(eigen(2)-lambda)*(eigen(3)-lambda))
   !
   !!write(*,*) "ASPH", delta, esse, rog
   !!!
   ! Store the info to get the average \Delta S and rog for all the nuclei, within the bulk and at the surface 
   delta_AVE=delta_AVE+delta 
   esse_AVE=esse_AVE+esse
   rog_AVE=rog_AVE+rog
   !
   !
   !!! JPCL 
   !!open(unit=5251, file='jpcl_alpha_1', status='unknown', position="append")
   !!open(unit=5252, file='jpcl_alpha_2', status='unknown', position="append")
   !!open(unit=5253, file='jpcl_alpha_3', status='unknown', position="append")
   !!open(unit=5254, file='jpcl_alpha_4', status='unknown', position="append")
   !!!
   !
   !
   if (surfF.eq.0) then ! within the bulk
      delta_AVE_BULK=delta_AVE_BULK+delta 
      esse_AVE_BULK=esse_AVE_BULK+esse
      rog_AVE_BULK=rog_AVE_BULK+rog
      !! JPCL stuff
      !if (mflag(counter+1).eq.1) then ! this one made it through
      !   write(5251,*) counter, delta
      !endif !else
      !   write(5252,*) counter, delta
      !!endif
      !!
   else ! at the surface
      delta_AVE_SURF=delta_AVE_SURF+delta
      esse_AVE_SURF=esse_AVE_SURF+esse
      rog_AVE_SURF=rog_AVE_SURF+rog
      !! JPCL stuff
      !if (mflag(counter+1).eq.1) then ! this one made it through
      !   write(5253,*) counter, delta
      !endif !else
      !   write(5254,*) counter, delta
      !!endif
      !!
   endif
   
   !! JPCL
   !close(5251)
   !close(5252)
   !close(5253)
   !close(5254)
endif


! Write down the colors for VMD
write(104,"("//adjustl(natformat)//"i10)") (r_color(k), k=1,nat)

! Write down rings statistics
write(stat_format,*) maxr-2
write(107,"(1E10.4,"//adjustl(stat_format)//"i10)") time, stat_nr(:)
if (trim(adjustl(switch_hbck)).eq.'yes') then
   write(307,"(1E10.4,"//adjustl(stat_format)//"i10)") time, stat_nr_HB(:)
endif
stat_nr_AVE(:)=stat_nr_AVE(:)+real(stat_nr(:))
stat_nr_HB_AVE(:)=stat_nr_HB_AVE(:)+real(stat_nr_HB(:))
deallocate(stat_nr,stat_nr_HB)

return

end subroutine rings

recursive function explore(index) result(fat)

use dfs
  implicit none
  integer ,intent(in) :: index
  integer :: j,iat
  integer :: fat


  followgraph(index)=1
  do j=1,neigh(index)
     iat=graph_solid_connect(index,j)
     if(followgraph(iat).eq.0)then
        predecessor(iat)=index
        followgraph(iat)=explore(iat)
     end if
  end do
  black=black+1
  fat=2

  ! who's who?
  lwho(count_cls,black)=cr_list(index)

end function explore

subroutine sort2(dati, n) ! Insertion sort

  integer :: n
  integer, dimension(n) :: dati
  integer :: i, j, tmp, min, pos

  do i=1,n-1
     min=dati(i)
     pos=i
     do j=i+1,n
        if (dati(j)>min) then
           min=dati(j) 
           pos=j
        end if
     end do
     tmp=dati(i) 
     dati(i)=dati(pos)
     dati(pos)=tmp
  end do

end subroutine sort2

! Find partcages 5^3 and 5^2 6
subroutine clath_cages(stat_wr,stat_nr,time,nat,natformat,kto)
    
    use MOD_vector3
    implicit none
    
    type(ragged_array) :: stat_wr
    integer, allocatable :: stat_nr(:) ! No. of 3,4,5,6,... rings
    
    type(vector3), allocatable :: rings_555(:), rings_655(:)
    integer :: n_rings_555, n_rings_655, nat, clath_color(nat), i, j, k
    real :: time
    character*100 :: natformat
    
    integer, dimension(:,:), allocatable :: rings5, rings6
    integer :: nrings5, nrings6, tmp_ring
    integer, allocatable :: n_cnx_55(:,:), n_cnx_65(:,:), kto(:)
    integer, allocatable :: t_n_cnx_55(:), t_n_cnx_56(:), t_n_cnx_6(:)
    
    type(cnx_graph), allocatable :: ring_cnxs_55(:), ring_cnxs_65(:)
    
    ! NOTE: stat_wr is an array of integer 2D arrays (stat_wr_size)
    ! The array is indexed by #members in ring, from 3 to max_rings (3->5, 4->6)
    ! 2D arrays are (which ring, member of ring), [#N-membered rings] x [N]
    
    rings5 = stat_wr%stat_wr_size(3)%mrings
    rings6 = stat_wr%stat_wr_size(4)%mrings
    nrings5 = stat_nr(3)
    nrings6 = stat_nr(4)
    
    call ringpairs(rings5,rings6,nrings5,nrings6,ring_cnxs_55,ring_cnxs_65,n_cnx_55,n_cnx_65, &
                   t_n_cnx_55,t_n_cnx_56,t_n_cnx_6)
    
    call partcage555(rings5,nrings5,ring_cnxs_55,n_cnx_55,t_n_cnx_55,n_rings_555,rings_555)
    
    call partcage655(rings5,rings6,nrings5,nrings6,ring_cnxs_55,ring_cnxs_65,n_cnx_55,n_cnx_65, &
                     t_n_cnx_55,t_n_cnx_56,t_n_cnx_6,n_rings_655,rings_655)
    
    deallocate(n_cnx_55, n_cnx_65, t_n_cnx_55, t_n_cnx_56, t_n_cnx_6, ring_cnxs_55, ring_cnxs_65)
    
    write(210,'(1E12.6,4X,I10,12X,I10)') time, n_rings_555, n_rings_655
    
    clath_color(:) = 0
    do i=1,n_rings_555 ; do j=1,3
        tmp_ring = rings_555(i)%rings(j)
        do k=1,5
            clath_color(kto(rings5(tmp_ring,k))) = 1
            clath_color(kto(rings5(tmp_ring,k))) = 1
        end do
    end do ; end do
    do i=1,n_rings_655
        tmp_ring = rings_655(i)%rings(1)
        do k=1,6
                if (clath_color(kto(rings6(tmp_ring,k))).eq.0) then
                    clath_color(kto(rings6(tmp_ring,k))) = 2
                else if (clath_color(kto(rings6(tmp_ring,k))).eq.1) then
                    clath_color(kto(rings6(tmp_ring,k))) = 3
                end if
        end do        
        do j=2,3
            tmp_ring = rings_655(i)%rings(j)
            do k=1,5
                if (clath_color(kto(rings5(tmp_ring,k))).eq.0) then
                    clath_color(kto(rings5(tmp_ring,k))) = 2
                else if (clath_color(kto(rings5(tmp_ring,k))).eq.1) then
                    clath_color(kto(rings5(tmp_ring,k))) = 3
                end if
            end do
        end do
    end do
    write(211,'('//adjustl(natformat)//'(I1,X))') (clath_color(i), i=1,nat)

end subroutine clath_cages

! Find 55 and 65 ring pairs
subroutine ringpairs(rings5,rings6,nrings5,nrings6,ring_cnxs_55,ring_cnxs_65,n_cnx_55,n_cnx_65, &
                     t_n_cnx_55,t_n_cnx_56,t_n_cnx_6)
    
    use MOD_vector3
    implicit none
    
    integer, dimension(:,:), allocatable :: rings5, rings6
    integer :: nrings5, nrings6, r1, r2, r3, o1, o2, i, j
    integer, allocatable :: n_cnx_55(:,:), n_cnx_65(:,:)
    integer, allocatable :: t_n_cnx_55(:), t_n_cnx_56(:), t_n_cnx_6(:)
    type(cnx_graph), allocatable :: ring_cnxs_55(:), ring_cnxs_65(:)
    
    allocate(n_cnx_55(nrings5,nrings5), n_cnx_65(nrings6,nrings5))
    allocate(t_n_cnx_55(nrings5), t_n_cnx_56(nrings5), t_n_cnx_6(nrings6))
    allocate(ring_cnxs_55(nrings5), ring_cnxs_65(nrings6))
    do i=1,nrings5
        allocate(ring_cnxs_55(i)%ring_cnx(nrings5))
    end do
    do i=1,nrings6
        allocate(ring_cnxs_65(i)%ring_cnx(nrings5))
    end do
    
    n_cnx_55(:,:) = 0
    n_cnx_65(:,:) = 0
    t_n_cnx_55(:) = 0
    t_n_cnx_56(:) = 0
    t_n_cnx_6(:) = 0
    
    do r1=1,nrings5
        if (r1.ne.nrings5) then ; do r2=r1,nrings5 ; do o1=1,5 ; do o2=1,5
            if (rings5(r1,o1).eq.rings5(r2,o2)) then
                n_cnx_55(r1,r2) = n_cnx_55(r1,r2) + 1
                t_n_cnx_55(r1) = t_n_cnx_55(r1) + 1
                t_n_cnx_55(r2) = t_n_cnx_55(r2) + 1
                ring_cnxs_55(r1)%ring_cnx(r2)%matches(n_cnx_55(r1,r2))%atom_match = (/ o1, o2 /)
            end if
        end do ; end do ; end do ; end if
        
        do r2=1,nrings6; do o1=1,5 ; do o2=1,6
            if (rings5(r1,o1).eq.rings6(r2,o2)) then
                n_cnx_65(r2,r1) = n_cnx_65(r2,r1) + 1
                t_n_cnx_56(r1) = t_n_cnx_56(r1) + 1
                t_n_cnx_6(r2) = t_n_cnx_6(r2) + 1
                ring_cnxs_65(r2)%ring_cnx(r1)%matches(n_cnx_65(r2,r1))%atom_match = (/ o2, o1 /)
            end if
        end do ; end do ; end do
    end do
    
end subroutine ringpairs

! Find partcages 555
subroutine partcage555(rings5,nrings5,ring_cnxs_55,n_cnx_55,t_n_cnx_55,n_rings_555,rings_555)
    
    use MOD_vector3
    implicit none
    
    integer, dimension(:,:), allocatable :: rings5
    integer :: nrings5, r1, r2, r3, o1, o2, i, j
    integer, allocatable :: n_cnx_55(:,:), t_n_cnx_55(:)
    type(cnx_graph), allocatable :: ring_cnxs_55(:)
    
    type(vector3), allocatable :: rings_555(:)
    integer :: n_rings_555
    
    allocate(rings_555(nrings5*(nrings5-1)*(nrings5-2)/6))
    
    ! We are looking for three rings with one common element,
    ! and one additional common element between each pair of rings
    
    n_rings_555 = 0
    
    do r1=1,nrings5-2
        ! First, check whether r1 is a possible candidate (i.e. has at least 4 connections)
        if (t_n_cnx_55(r1).ge.4) then ; do r2=r1+1,nrings5-1
            if ((t_n_cnx_55(r2).ge.4).and.(n_cnx_55(r1,r2).eq.2)) then ; do r3=r2+1,nrings5
                if ((n_cnx_55(r1,r3).eq.2).and.(n_cnx_55(r2,r3).eq.2)) then
                    ! Now have three rings with two connections each. Must check one is common.
                    outer: do i=1,2 ; do j=1,2
                        if (ring_cnxs_55(r1)%ring_cnx(r2)%matches(i)%atom_match(1).eq.&
                            &ring_cnxs_55(r1)%ring_cnx(r3)%matches(j)%atom_match(1)) then
                            n_rings_555 = n_rings_555 + 1
                            rings_555(n_rings_555)%rings = (/ r1, r2, r3 /)
                            exit outer
                        end if
                    end do ; end do outer
                end if
            end do ; end if
        end do ; end if
    end do

end subroutine partcage555


! Find partcages 655
subroutine partcage655(rings5,rings6,nrings5,nrings6,ring_cnxs_55,ring_cnxs_65,n_cnx_55,n_cnx_65, &
                       t_n_cnx_55,t_n_cnx_56,t_n_cnx_6,n_rings_655,rings_655)
    
    use MOD_vector3
    implicit none
    
    integer, dimension(:,:), allocatable :: rings5, rings6
    integer :: nrings5, nrings6, r1, r2, r3, o1, o2, i, j
    integer, allocatable :: n_cnx_55(:,:), n_cnx_65(:,:)
    integer, allocatable :: t_n_cnx_55(:), t_n_cnx_56(:), t_n_cnx_6(:)
    type(cnx_graph), allocatable :: ring_cnxs_55(:), ring_cnxs_65(:)
    
    type(vector3), allocatable :: rings_655(:)
    integer :: n_rings_655
    
    allocate(rings_655(nrings6*nrings5*(nrings5-1)/2))
    
    ! We are looking for three rings (2x5 & 1x6) with one common element,
    ! and one additional common element between each pair of rings
    
    n_rings_655 = 0
    
    do r1=1,nrings6
        ! First, check whether r1 is a possible candidate (i.e. has at least 4 connections)
        if (t_n_cnx_6(r1).ge.4) then ; do r2=1,nrings5-1
            if ((t_n_cnx_55(r2).ge.2).and.(t_n_cnx_56(r2).ge.2).and.(n_cnx_65(r1,r2).eq.2)) then ; do r3=r2+1,nrings5
                if ((n_cnx_65(r1,r3).eq.2).and.(n_cnx_55(r2,r3).eq.2)) then
                    ! Now have three rings with two connections each. Must check one is common.
                    outer: do i=1,2 ; do j=1,2
                        if (ring_cnxs_65(r1)%ring_cnx(r2)%matches(i)%atom_match(1).eq.&
                            &ring_cnxs_65(r1)%ring_cnx(r3)%matches(j)%atom_match(1)) then
                            n_rings_655 = n_rings_655 + 1
                            rings_655(n_rings_655)%rings = (/ r1, r2, r3 /)
                            exit outer
                        end if
                    end do ; end do outer
                endif
            end do ; end if
        end do ; end if
    end do

end subroutine partcage655


end module MOD_rings
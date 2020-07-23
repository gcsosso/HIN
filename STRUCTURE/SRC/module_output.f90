module MOD_output

contains

subroutine output(dostuff,lframe,fframe,stride,switch_outxtc,ns,ws,n_ws,zmesh,dens,nz,dz,box_trans, &
                  switch_rings,r_ns,r_ws,n_r_ws,maxr,stat_nr_AVE,switch_cages,n_ddc_AVE,n_hc_AVE, &
                  switch_hex,n_hex_AVE,switch_bonds,npairs,nz_bAVE,b_zmin,b_dz,b_bins,b_bmax,b_bmin, &
                  pdbon_AVE,npairs_cn,cn_AVE,switch_cls,n_cls_AVE,cart,switch_zdens, &
                  switch_r_cls,r_cls_W,switch_xyfes,xydens,xymax,xymin,nxy,xmesh,ymesh,nsurf,nbulk,n_ddc_AVE_SURF, &
                  n_hc_AVE_SURF,n_hex_AVE_SURF, &
                  n_ddc_AVE_BULK,n_hc_AVE_BULK,n_hex_AVE_BULK,switch_ffss, &
                  delta_AVE,delta_AVE_BULK,delta_AVE_SURF,esse_AVE,esse_AVE_BULK,esse_AVE_SURF, &
                  rog_AVE,rog_AVE_BULK,rog_AVE_SURF,ze_AVE,ze_AVE_BULK, &
                  ze_AVE_SURF,d_charge,switch_electro,e_nz,e_zmesh, &
                  switch_th,switch_water,o_nz,o_zmesh,w_order,zop_AVE,stat_nr_HB_AVE,switch_hbck)

implicit none

! Local
integer :: i, j, k, ibin, l, n
real :: rstep, h, rsum
real, parameter :: epsi=0.0055267840353714 ! permettivity of vacuum in e/(V*angs)
real, allocatable :: efield(:), epot(:)
character*100 :: wformat

! Arguments
integer :: dostuff, fframe, stride, lframe, nz, b_bins, nz_bAVE, e_nz, o_nz
integer :: ns, r_ns, npairs, npairs_cn, maxr, cart, nxy, nsurf, nbulk
integer, allocatable :: n_ws(:), n_r_ws(:)
real :: box_trans(cart,cart), zmin, zmax, r_zmin, r_zmax, dz, zop_AVE
real :: b_zmin, b_zmax, b_dz, b_bmin, b_bmax, xymax, xymin, ze_AVE, ze_AVE_BULK, ze_AVE_SURF
real :: n_ddc_AVE_SURF, n_hc_AVE_SURF, n_hex_AVE_SURF, n_ddc_AVE_BULK, n_hc_AVE_BULK, n_hex_AVE_BULK
real :: n_ddc_AVE, n_hc_AVE, n_hex_AVE, n_cls_AVE
real :: delta_AVE, delta_AVE_BULK, delta_AVE_SURF, esse_AVE, esse_AVE_BULK, esse_AVE_SURF, rog_AVE, rog_AVE_BULK, rog_AVE_SURF
real, allocatable :: dens(:,:), zmesh(:), stat_nr_AVE(:), pdbon_AVE(:,:,:), cn_AVE(:,:), stat_nr_HB_AVE(:)
real, allocatable :: xydens(:,:,:), xmesh(:), ymesh(:), d_charge(:), e_zmesh(:), o_zmesh(:), w_order(:)
logical(1) :: switch_outxtc, switch_zdens, switch_rings, switch_cls, switch_bonds, switch_xyfes, switch_hbck
logical(1) :: switch_hex, switch_cages, switch_r_cls, switch_ffss, switch_electro, switch_th
character(3) :: switch_water, r_cls_W
character(4), allocatable :: ws(:)
character(5), allocatable :: r_ws(:)

if (dostuff.ne.((lframe-fframe)/stride)+1) then
   write(99,*) "Something's wrong with the number of frames!"
   write(99,*) "We've taken into account ", dostuff, "frames, but you have specified ", ((lframe-fframe)/stride)+1
   stop
else
   write(99,*) "We have taken into account ", dostuff, "frames, from frame ", fframe, " to frame ", lframe, " with stride ", stride
endif

if (switch_outxtc) then
   write(99,*) "We have written a subset of the original .xtc trajectory. See: hin_structure.out.xtc"
endif

if (switch_zdens) then
   write(99,*) "We have computed the number density profile along z. See: hin_structure.out.zdens"      
   write(99,*) "We have ", ns, "atomic species, namely ", (ws(i), i=1,ns)
   do j=1,ns
      write(99,*) "We have ", n_ws(j), ws(j)
   enddo
   open(unit=102, file='hin_structure.out.zdens', status='unknown')
   ! Output the number density [n. atoms / angs^3] as a function of z [angs]
   write(wformat,*) ns+1
   zmesh(:)=zmesh(:)*10.0d0
   dens(:,:)=(dens(:,:)/real(dostuff))/(box_trans(1,1)*box_trans(2,2)*dz*1000.0)
   do k=1,nz
      write(102,"("//adjustl(wformat)//"f20.10)") zmesh(k), dens(k,:)
   enddo
endif

if (switch_xyfes) then
   write(99,*) "We have computed the 2D FES in the xy plane. See: hin_structure.out.xyfes"
   write(99,*) "We have ", ns, "atomic species, namely ", (ws(i), i=1,ns)
   do j=1,ns
      write(99,*) "We have ", n_ws(j), ws(j)
   enddo
   open(unit=1102, file='hin_structure.out.xyfes', status='unknown')
   ! Output the 2D fes [n. atoms / angs^3] as a function of x and y [angs]
   write(wformat,*) ns+2
   xmesh(:)=xmesh(:)*10.0
   ymesh(:)=ymesh(:)*10.0
   xydens(:,:,:)=(xydens(:,:,:)/real(dostuff))/(box_trans(1,1)*box_trans(2,2)*(abs(xymax-xymin))*1000.0)
   do k=1,nxy
      write(1102,*)
      do l=1,nxy
         write(1102,*) xmesh(k), ymesh(l),  xydens(k,l,:) 
      enddo
   enddo
endif

if (switch_rings) then
   write(99,*) "We have computed rings statistics. See: hin_structure.out.rings.stats"
   write(99,*) "We have ", r_ns, "atomic species, namely ", (r_ws(i), i=1,r_ns)
   do j=1,r_ns
      write(99,*) "We have ", n_r_ws(j), r_ws(j)
   enddo
   do n=3,maxr
      write(99,*) "Average number of ", n , "-membered rings = ", stat_nr_AVE(n)/real(dostuff)
   enddo
   if (switch_hbck) then
      write(99,*) "We have also computed the number or wholly hydrogen-bonded rings!!"
      write(99,*) "See: hin_structure.out.rings.stats.HB"
      write(99,*) "Here are the averages ** for the HB rings only **"
      do n=3,maxr
         write(99,*) "Average number of ", n , "-membered rings = ", stat_nr_HB_AVE(n)/real(dostuff)
      enddo
   endif

   if (switch_cages) then
      write(99,*) "Average number of DDCcages: ", n_ddc_AVE/real(dostuff)
      write(99,*) "Average number of HCcages: ", n_hc_AVE/real(dostuff)
   endif
   if (switch_ffss) then 
      if (nbulk+nsurf.ne.dostuff) write(99,*) "Something's wrong, mate!"
      write(99,*) "Of the ", nbulk+nsurf, " clusters at this interface, "
      write(99,*) nbulk, " are within the bulk and "
      write(99,*) nsurf, " are on top of the surface"
      write(99,*) "Average number of DDCcages within the BULK: ", n_ddc_AVE_BULK/real(nbulk)
      write(99,*) "Average number of HCcages within the BULK: ", n_hc_AVE_BULK/real(nbulk)
      write(99,*) "Average number of DDCcages within the SURFACE: ", n_ddc_AVE_SURF/real(nsurf)
      write(99,*) "Average number of HCcages within the SURFACE: ", n_hc_AVE_SURF/real(nsurf)
      write(99,*) "Asphericity"
      write(99,*) "Average on all the clusters (delta, esse, rog): ", &
                  delta_AVE/real(dostuff), esse_AVE/real(dostuff) , rog_AVE/real(dostuff)
      write(99,*) "Asph BULK: ", delta_AVE_BULK/real(nbulk), esse_AVE_BULK/real(nbulk) , rog_AVE_BULK/real(nbulk)
      write(99,*) "Asph SURF: ", delta_AVE_SURF/real(nsurf), esse_AVE_SURF/real(nsurf) , rog_AVE_SURF/real(nsurf)
      write(99,*) "Z_EXT all: ", ze_AVE/real(dostuff)
      write(99,*) "Z_EXT BULK: ", ze_AVE_BULK/real(nbulk)
      write(99,*) "Z_EXT SURF: ", ze_AVE_SURF/real(nsurf)
   endif
   if (switch_hex) then
      write(99,*) "Average number of hexagonal rings: ", n_hex_AVE/real(dostuff)
   endif
   ! Cluster hexagonal rings, e.g. to find the largest patch of hexagonal rings sitting on top of the surface
   if (switch_r_cls) then
      if (trim(adjustl(r_cls_W)).ne.'SIX') then
         write(99,*) "Sorry mate, I can do only six membered rings at the moment..."
      else
      write(99,*) "We have also clustered hexagonal rings together..."
      write(99,*) "See hin_structure.out.rings.stats.patch for the n. of atoms involved in the biggest hexagonal patch..."
      endif
   endif
endif

if (switch_bonds) then
   write(99,*) "We have computed the probability density distribution of bond lengths. See: hin_structure.out.bonds"
   write(99,*) "We have ", ns, "atomic species, namely ", (ws(i), i=1,ns)
   write(99,*) "We have ", npairs, " bond pairs"
   open(unit=105, file='hin_structure.out.bonds', status='unknown')
   ! Output the probability density distribution of bond lengths as a function of z [angs]
   do k=1,nz_bAVE
      zmesh(k)=(b_zmin+((k*b_dz)-(b_dz/2.0)))*10.0
   enddo
   write(wformat,*) npairs+2
   do k=1,nz_bAVE
      write(105,*) " "
      write(105,*) " "
      do ibin=1,b_bins
         rstep=(((abs(b_bmax-b_bmin)/real(b_bins))*dble(ibin))+b_bmin)*10.0
         write(105,"("//adjustl(wformat)//"f20.10)") rstep, (pdbon_AVE(k,i,ibin)/real(dostuff), i=1,npairs), zmesh(k)
      enddo
   enddo
   write(99,*) "We have computed the coordination number. See: hin_structure.out.coord"
   write(99,*) "We have ", ns, "atomic species, namely ", (ws(i), i=1,ns)
   write(99,*) "We have ", npairs_cn, " cn pairs"
   ! which ones? write them down !
   open(unit=106, file='hin_structure.out.coord', status='unknown')
   ! Output the coordination number
   write(wformat,*) npairs_cn+1
   do k=1,nz_bAVE
      write(106,"("//adjustl(wformat)//"f20.10)") (cn_AVE(k,i)/real(dostuff), i=1,npairs_cn), zmesh(k)
   enddo
endif

if (switch_cls) then
   write(99,*) "We got the biggest ice-like cluster in time. See: hin_structure.out.cls.lambda"
   write(99,*) "Average nucleus size: ", n_cls_AVE/real(dostuff)
   write(99,*) "Colors are in hin_structure.out.cls.color"
endif

if (switch_electro) then
   write(99,*) "We have calculated some electrostatic as well. See: hin_structure.out.electro"
   open(unit=902, file='hin_structure.out.electro', status='unknown')

   ! Output:
   ! 1. z [angs] 
   ! 2. charge density profile  [e/angs^3]
   ! 3. electric field profile  [V/angs]
   ! 4. electrostatic potential [V]

   e_zmesh(:)=e_zmesh(:)*10.0
   d_charge(:)=(d_charge(:)/real(dostuff))/(box_trans(1,1)*box_trans(2,2)*dz*1000.0)

   ! Integrate the charge to get the electric field
   allocate(efield(e_nz))
   efield(:)=0.0
   rsum=0.0
   h=((e_zmesh(2)-e_zmesh(1)))/(2.0*epsi)
   do k=1,e_nz-1
      rsum=rsum+((d_charge(k)+d_charge(k+1))*h)
      efield(k)=rsum
   enddo 

   ! Integrate the field to get the electrostatic potential
   allocate(epot(e_nz))  
   epot(:)=0.0
   rsum=0.0
   h=((e_zmesh(2)-e_zmesh(1)))
   do k=1,e_nz-2
      rsum=rsum-((efield(k)+efield(k+1))*h)
      epot(k)=rsum
   enddo

   do k=1,e_nz-2
      write(902,'(4f20.10)') e_zmesh(k), d_charge(k), efield(k), epot(k)
   enddo
   close(902)
endif

if (switch_th) then
   write(99,*) "We have calculated some order parameters as well..."
   if (trim(adjustl(switch_water)).eq.'yes') then
      write(99,*) "Water ordering - profile along z..."
      open(unit=254, file='hin_structure.out.w_order', status='unknown')
      o_zmesh(:)=o_zmesh(:)*10.0
      w_order(:)=(w_order(:)/real(dostuff))!/(box_trans(1,1)*box_trans(2,2)*dz*1000.0)      
      do k=1,o_nz
         write(254,'(2f20.10)') o_zmesh(k), w_order(k)
      enddo
   else
      write(99,*) "Average angle of the molecular axis of choice wrt the z-axis"
      write(99,*) zop_AVE/real(dostuff)
   endif

endif

end subroutine output

end module MOD_output

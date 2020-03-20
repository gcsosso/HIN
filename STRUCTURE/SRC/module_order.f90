module MOD_order

contains

subroutine order_alloc(o_nz,o_zmax,o_zmin,o_dz,w_order,o_zmesh,switch_water)

implicit none

! Arguments

integer :: o_nz
real :: o_zmax, o_zmin, o_dz, middle
real, allocatable :: w_order(:), o_zmesh(:)
character*3 :: switch_water

! Local

integer :: k

! Allocate the d_charge array...
o_nz=int((o_zmax-o_zmin)/o_dz)
if (o_nz.lt.1) then
   write(99,*) "Your mesh along z (order section) is ill-defined. Modify o_zmin, o_zmax and/or o_dz"
   stop
endif
allocate(w_order(o_nz),o_zmesh(o_nz))
w_order(:)=0.0
do k=1,o_nz
   o_zmesh(k)=o_zmin+((k*o_dz)-(o_dz/2.0))
enddo

if (trim(adjustl(switch_water)).eq.'mol') then
	open(unit=254, file='hin_structure.out.w_order', status='unknown')
	open(unit=255, file='hin_structure.out.w_order.color', status='unknown')
endif

end subroutine order_alloc

subroutine order(o_nz,o_zmax,o_zmin,o_dz,w_order,o_zmesh,nat,pos, &
                mq_all,cart,middle,switch_water,sym,wmol,resname, &
                resnum,axis_1,axis_2,zop_AVE,natformat,icell)

implicit none

!! Arguments
integer :: o_nz, cart,nat
real :: o_zmax, o_zmin, o_dz, middle
real, parameter :: rad2deg=57.2958
real, allocatable :: w_order(:), o_zmesh(:)
real, allocatable :: mq_all(:), pos(:,:)
character*3 :: switch_water
character*4 :: wmol, axis_1, axis_2
character*4, allocatable :: sym(:)
character*5, allocatable :: resname(:)
integer, allocatable :: resnum(:)
character*100 :: natformat
real :: icell(cart*cart)

! Local
integer :: i, j, k, flag_1, flag_2, idx_1, idx_2, resn_1, resn_2
real, allocatable :: nslice(:), wok(:)
real :: lb, ub, com(cart), mmm, dm(cart), pos_1(cart), pos_2(cart), zop, nmol, zop_AVE
character*4 :: sym_1, sym_2
real, allocatable :: w_order_col(:), w_order_mol(:), w_oz(:)
integer :: n_mol
character*100 :: n_mol_format

! Remove COM
com(:)=0.0
mmm=0.0
do j=1,nat
   com(cart)=com(cart)+(pos(cart,j)*mq_all(j))
   mmm=mmm+mq_all(j)
enddo
com(cart)=com(cart)/mmm

if (trim(adjustl(switch_water)).eq.'yes') then
 
   allocate(nslice(o_nz),wok(o_nz)) ; nslice(:)=0.0 ; wok(:)=0.0
   !! Orientational order parameter (wrt to the normal to the slab)
   !! as a function of z (water only!)
   do j=1,nat
      ! rescale the positions wrt the center of mass of each frame - than move the whole thing in the middle
      pos(cart,j)=pos(cart,j)-com(cart)+middle
      do k=1,o_nz
         lb=o_zmin+(k-1)*o_dz
         ub=lb+o_dz
         ! save time by ignoring non-water stuff (and HW and MW as well!)
         if (trim(adjustl(sym(j))).ne.trim(adjustl(axis_1))) cycle 
         ! build the profile according to the position of the oxygens...
         if (pos(cart,j).gt.lb.and.pos(cart,j).le.ub) then
         nslice(k)=nslice(k)+1.0 ! number of molecules within the slice...
         ! get the unit vector along the dipole moment
         ! which we assume it lies along the TIP4P OW-MW segment (bisector of the H-O-H angle)
         ! there should be no need of invoking pbc - even if the molecules are not whole
         dm(:)=pos(:,j)+pos(:,j)-pos(:,j+1)-pos(:,j+2) ! MW ----> OW - chemistry convention (physicist would do the other way around
                                   ! for the dipole moment, e.g. - -> + instead of + -> - )
         dm(:)=dm(:)/(sqrt(dm(1)**2.0+dm(2)**2.0+dm(3)**2.0)) ! From -1 to 1
         !(acos(dm(cart)))*rad2deg ! angle between the dipole moment of the water molecule and the z-axis
         ! Ranges from 0 deg (dm(cart=1)) to 180 deg (dm(cart)=-1). Average is 90 deg. 
         ! order parameter as a function of z. average over all the molecules within a slice
         ! angle
         wok(k)=wok(k)+(acos(dm(cart)))*rad2deg
         ! Values > 90 deg correspond to negative values of dm(cart).
         ! That is, negative value of the projection of the (water) dipole moment on the z-axis
         ! That is, the water dipole moment is pointing down - with respect to z
         ! Values < 90 - the molecular axis or dipole moment or whatever is pointing in the same dir as z
         endif
      enddo
   enddo

   ! normalize by the number of molecules in each slice
   do k=1,o_nz
      !if (nslice(k).gt.0.0) then
      !   wok(k)=wok(k)/nslice(k)
      !else
      !   wok(k)=0.0
      !endif
      wok(k)=wok(k)/nslice(k)
      !if (nslice(k).gt.0.0) then
      !   write(*,*) wok(k)
      !endif
      w_order(k)=w_order(k)+wok(k) ! pile up, to be averaged in the output...
   enddo

   ! Get also the Probability density of these angles for all the occurrences of 
   ! water molecules within the slabs closest to e.g. CHL ! 

else if (trim(adjustl(switch_water)).eq.'mol') then ! We are outputing an order parameter for each water molecule

   allocate(w_order_col(nat), w_order_mol(nat), w_oz(nat))
   w_order_col(:) = 0.0
   w_order_mol(:) = 0.0
   w_oz(:) = 0.0
   n_mol = 0
	
   !! Orientational order parameter (wrt to the normal to the slab)
   !! to be written to a color file (water only!)
   do j=1,nat
	 	! save time by ignoring non-water stuff (and HW and MW as well!)
	 	if (trim(adjustl(sym(j))).ne.trim(adjustl(axis_1))) cycle
		if (pos(cart,j).gt.o_zmin.and.pos(cart,j).le.o_zmax) then ! filter to just the desired slice
            n_mol = n_mol+1
            w_oz(n_mol) = pos(cart,j)
	 	    ! get the unit vector along the dipole moment
	 	    ! which we assume it lies along the TIP4P OW-MW segment (bisector of the H-O-H angle)
	 	    ! there should be no need of invoking pbc - even if the molecules are not whole
            dm(:)=pos(:,j)+pos(:,j)-pos(:,j+1)-pos(:,j+2) ! MW ----> OW - chemistry convention (physicist would do the other way around
            call images(cart,0,1,1,icell,dm(1),dm(2),dm(3))
										  ! for the dipole moment, e.g. - -> + instead of + -> - )
	 	    dm(:)=dm(:)/(sqrt(dm(1)**2.0+dm(2)**2.0+dm(3)**2.0)) ! From -1 to 1
	 	    !(acos(dm(cart)))*rad2deg ! angle between the dipole moment of the water molecule and the z-axis
		    ! Ranges from 0 deg (dm(cart=1)) to 180 deg (dm(cart)=-1). Average is 90 deg. 
	 	    ! order parameter as a function of z. average over all the molecules within a slice
	 	    ! angle
	 	    w_order_col(j) = (acos(dm(cart)))*rad2deg
	 	    w_order_mol(n_mol) = (acos(dm(cart)))*rad2deg
        end if
	 	! Values > 90 deg correspond to negative values of dm(cart).
	 	! That is, negative value of the projection of the (water) dipole moment on the z-axis
	 	! That is, the water dipole moment is pointing down - with respect to z
	 	! Values < 90 - the molecular axis or dipole moment or whatever is pointing in the same dir as z
   enddo
	write(n_mol_format,*) n_mol
	
	write(255,'('//adjustl(natformat)//'F11.4)') (w_order_col(i), i=1,nat)
	
	write(254,'('//adjustl(n_mol_format)//'F11.4)') (w_oz(i), i=1,n_mol)
	write(254,'('//adjustl(n_mol_format)//'F11.4)') (w_order_mol(i), i=1,n_mol)
	

else ! molecules other than water

   !! Orientational order parameter (wrt to the normal to the slab)
   flag_1=0 
   flag_2=0
   zop=0.0
   nmol=0.0

   ! Distribution - to be fixed TBF !
   !open(unit=845, file='oh_d.dat', status='unknown')
   ! END TBF
 
   do j=1,nat
      ! rescale the positions wrt the center of mass of each frame - than move the whole thing in the middle
      pos(cart,j)=pos(cart,j)-com(cart)+middle
      if (trim(adjustl(sym(j))).eq.trim(adjustl(axis_1)).and.pos(cart,j).gt.o_zmin.and.pos(cart,j).lt.o_zmax) then
         pos_1(:)=pos(:,j)
         flag_1=1 ; resn_1=resnum(j); idx_1=j; sym_1=trim(adjustl(sym(j)))
      endif
      if (trim(adjustl(sym(j))).eq.trim(adjustl(axis_2)).and.pos(cart,j).gt.o_zmin.and.pos(cart,j).lt.o_zmax) then
         pos_2(:)=pos(:,j)
         flag_2=1 ; resn_2=resnum(j); idx_2=j; sym_2=trim(adjustl(sym(j))) 
      endif
   !   ! All of them
      if (flag_1.eq.1.and.flag_2.eq.1.and.resn_1.eq.resn_2) then
   !   ! Silly hardcoded one to consider some z-region...TBF !
   !    if (flag_1.eq.1.and.flag_2.eq.1.and.pos_1(cart).lt.10.0.and.pos_2(cart).lt.10.0) then ! Gets only the bottom layer...
      !

         !write(*,*) resn_1, resn_2, idx_1, idx_2, sym_1, sym_2, pos_1(3), pos_2(3)
     
         nmol=nmol+1.0
         dm(:)=pos_2(:)-pos_1(:) 
         !!! DEBUG
         !write(*,*) idx_1, sym_1, resn_1
         !write(*,*) idx_2, sym_2, resn_2
         !!! END DEBUG
         !
         dm(:)=dm(:)/(sqrt(dm(1)**2.0+dm(2)**2.0+dm(3)**2.0))
! why abs?!         zop=zop+(acos(abs(dm(cart))))*rad2deg ! angle between the molecular axis of choice and the z-axis
         ! Distribution - TBF
         !write(*,*) (acos(dm(cart)))*rad2deg
         ! END TBF
         ! up and down !
         zop=zop+(acos(dm(cart)))*rad2deg ! angle between the molecular axis of choice and the z-axis
         flag_1=0 ; flag_2=0
      endif

   enddo

   zop=zop/nmol
   zop_AVE=zop_AVE+zop
   !write(*,*) zop, zop_AVE, nmol

endif

end subroutine order

end module MOD_order

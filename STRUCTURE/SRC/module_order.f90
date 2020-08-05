module MOD_order

contains

subroutine order_alloc()

implicit none

 open(unit=254, file='hin_structure.out.w_order', status='unknown')
 open(unit=255, file='hin_structure.out.w_order.color', status='unknown')

end subroutine order_alloc

subroutine order(nat, pos, cart, icell, n_filtered, list_filtered, filt_param, switch_filt_param, natformat)

implicit none

!! Arguments
integer :: cart, nat, n_filtered, list_filtered(:)
real, parameter :: rad2deg=57.2958
real, allocatable :: w_order(:), o_zmesh(:)
real :: pos(:,:), filt_param(:)
character*100 :: natformat
real :: icell(cart*cart)
logical(1) :: switch_filt_param

! Local
integer :: i, ii
real :: dm(cart)
real :: w_order_col(nat), w_order_mol(n_filtered)
character*100 :: n_mol_format

w_order_col(:) = 0.0
w_order_mol(:) = 0.0

!! Orientational order parameter (wrt to the normal to the slab)
! Values > 90 deg correspond to negative values of dm(cart).
! That is, negative value of the projection of the (water) dipole moment on the z-axis
! That is, the water dipole moment is pointing down - with respect to z
! Values < 90 - the molecular axis or dipole moment or whatever is pointing in the same dir as z
do ii=1,n_filtered ; i = list_filtered(ii)
   ! get the unit vector along the dipole moment
   ! which we assume it lies along the TIP4P OW-MW segment (bisector of the H-O-H angle)
   ! there should be no need of invoking pbc - even if the molecules are not whole
   dm(:)=pos(:,i)+pos(:,i)-pos(:,i+1)-pos(:,i+2) ! MW ----> OW - chemistry convention (physicist would do the other way around
   call images(cart,0,1,1,icell,dm(1),dm(2),dm(3))
                                 ! for the dipole moment, e.g. - -> + instead of + -> - )
   dm(:)=dm(:)/(sqrt(dm(1)**2.0+dm(2)**2.0+dm(3)**2.0)) ! From -1 to 1
   !(acos(dm(cart)))*rad2deg ! angle between the dipole moment of the water molecule and the z-axis
   ! Ranges from 0 deg (dm(cart=1)) to 180 deg (dm(cart)=-1). Average is 90 deg. 
   ! order parameter as a function of z. average over all the molecules within a slice
   ! angle
   w_order_col(i) = (acos(dm(cart)))*rad2deg
   w_order_mol(ii) = (acos(dm(cart)))*rad2deg
enddo
 write(n_mol_format,*) n_filtered

 write(255,'('//adjustl(natformat)//'F11.4)') (w_order_col(i), i=1,nat)
 write(254,'('//adjustl(n_mol_format)//'I7)') (list_filtered(ii), ii=1,n_filtered)
 if (switch_filt_param) write(254,'('//adjustl(n_mol_format)//'F11.4)') (filt_param(ii), ii=1,n_filtered)
 write(254,'('//adjustl(n_mol_format)//'F11.4)') (w_order_mol(ii), ii=1,n_filtered)
	

end subroutine order

end module MOD_order

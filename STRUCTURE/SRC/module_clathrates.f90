module MOD_clathrates

contains

! Creates output directories and files for F3 and/or F4 calculations 
subroutine clathrates_alloc(switch_f3,switch_f4)

implicit none

character*3 :: switch_f3, switch_f4

! Make the tmp dir and open output files
call system("rm -r -f data-c ; mkdir data-c")
if (trim(adjustl(switch_f3)).eq.'yes') then
    write(99,*) "We are calculating the clathrate F3 order parameter. See hin_structure.out.clathrates.f3"
    open(unit=203, file='hin_structure.out.clathrates.f3', status='unknown')
    write(203,*) "Some sort of F3 order parameter statistics will be written here."
endif
if (trim(adjustl(switch_f4)).eq.'yes') then
    write(99,*) "We are calculating the clathrate F4 order parameter. See hin_structure.out.clathrates.f4"
    open(unit=204, file='hin_structure.out.clathrates.f3', status='unknown')
    write(204,*) "Some sort of F3 order parameter statistics will be written here."
endif
open(unit=205, file='hin_structure.out.clathrates.color', status='unknown')
open(unit=206, file='hin_structure.out.clathrates.stats', status='unknown')

end subroutine clathrates_alloc

subroutine clathrates_f3(f_zmin,f_zmax,f_cut,f_ns,f_ws,n_f_ws,list_f_ws, &
                        counter,kto,icell,box_trans)

implicit none

real :: f_zmin, f_zmax, f_cut, icell(cart*cart), box_trans(cart,cart)
integer :: f_ns, counter, nxyz
integer, allocatable :: n_f_ws(:), kto(:)
character*4, allocatable :: f_ws(:)

! Getting the - variable - box...
icell(1)=box_trans(1,1) ; icell(2)=box_trans(1,2) ; icell(3)=box_trans(1,3)
icell(4)=box_trans(2,1) ; icell(5)=box_trans(2,2) ; icell(6)=box_trans(2,3)
icell(7)=box_trans(3,1) ; icell(8)=box_trans(3,2) ; icell(9)=box_trans(3,3)

! Write down an .xyz with the region we have selected
open(unit=69, file='conf.xyz', status='unknown')
open(unit=70, file='tmp.dat', status='unknown')
nxyz=0
kto(:)=0
do i=1,f_ns
  do j=1,n_f_ws(i)
     if (pos(cart,list_f_ws(i,j)).ge.f_zmin.and.pos(cart,list_r_ws(i,j)).le.f_zmax) then
        nxyz=nxyz+1
        ! Index nxyz in conf.xyz corresponds to index list_r_ws(i,j) in the global .xtc
        ! Store this information for visualisation purposes
        kto(nxyz)=list_f_ws(i,j)
     endif
  enddo 
enddo 
write(70,*) nxyz
write(70,'(4f20.10)') icell(1)*10.0, icell(5)*10.0, icell(9)*10.0, f_cut*10.0
do i=1,f_ns
  do j=1,n_f_ws(i)
     if (pos(cart,list_f_ws(i,j)).ge.f_zmin.and.pos(cart,list_f_ws(i,j)).le.f_zmax) then
        !write(69,'(1a5,3f20.10)') sym(list_r_ws(i,j)), pos(:,list_r_ws(i,j))*10.0
        ! We write down everything as oxygen atoms!! 
        write(69,'(1a5,3f20.10)') "O", pos(:,list_f_ws(i,j))*10.0
     endif
  enddo
enddo
close(69)
close(70)


end subroutine clathrates_f3

end module MOD_rings
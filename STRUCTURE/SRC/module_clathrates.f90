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

subroutine clathrates_f3(f_zmin,f_zmax,f_cut,f_ns,f_ws,n_f_ws,list_f_ws,counter)

implicit none

real :: f_zmin, f_zmax, f_cut
integer :: f_ns, counter                ! Number of species of interest, frame
integer, allocatable :: n_f_ws(:)       ! Number of each species of interest (probably only OW)
integer, allocatable :: list_f_ws(:,:)  ! Atom numbers of species of interest
character*4, allocatable :: f_ws(:)     ! List of species of interest
integer :: first_coord_shell(10)        ! First coordination shell of the current atom
integer :: size_first_coord_shell       ! Size of first coordination shell
real :: dx, dy, dz                      ! X, Y and Z distances between two atoms

do i=1,f_ns ! Iterate through species of interest (probably just OW)
    do j=1,n_f_ws(i) ! Iterate through atoms of that species
        if (pos(cart,list_f_ws(i,j)).ge.f_zmin.and.pos(cart,list_f_ws(i,j)).le.f_zmax) then
            ! If atom is in Z-region of interest, calculate it's first coordination shell
            first_coord_shell(:) = 0
            size_first_coord_shell = 0
            do k=1,n_f_ws(i) ! Iterate through other atoms of species
                if (j/=k) then
                    dx = pos(1,list_f_ws(i,k))-pos(1,list_f_ws(i,j))
                    dy = pos(2,list_f_ws(i,k))-pos(2,list_f_ws(i,j))
                    dz = pos(3,list_f_ws(i,k))-pos(3,list_f_ws(i,j))
                    if (dx*dx + dy*dy + dz*dz <= f_cut*f_cut) then
                        size_first_coord_shell = size_first_coord_shell + 1
                        first_coord_shell(size_first_coord_shell) = list_f_ws(i,k)
                    endif
                endif
            enddo
        endif
    enddo
enddo
            

end subroutine clathrates_f3

end module MOD_rings
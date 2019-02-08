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

subroutine clathrates_f3(f_zmin,f_zmax,f_cut,f_ns,f_ws,n_f_ws,list_f_ws,counter, &
                         cart,icell)

implicit none

real :: f_zmin, f_zmax, f_cut
integer :: f_ns, counter                ! Number of species of interest, frame
integer, allocatable :: n_f_ws(:)       ! Number of each species of interest (probably only OW)
integer, allocatable :: list_f_ws(:,:)  ! Atom indices of species of interest
integer, allocatable :: tot_atoms(:)    ! Count of number of atoms for which F3 is calculated
character*4, allocatable :: f_ws(:)     ! List of species of interest
integer :: first_coord_shell(10,5)      ! First coordination shell of the current atom: (index, dx, dy, dz, dsq)
integer :: size_first_coord_shell       ! Size of first coordination shell
real :: dx, dy, dz                      ! X, Y and Z distances between two atoms
real :: dsq                             ! Square distance between two atoms
real :: F3_part, F3_atom                ! F3 parameter for triples, atoms
real, allocatable :: F3_avg(:)          ! F3 parameter for frame-wide avg
real :: kdotl, cos2_num, cos2_den       ! k.l, |cos|cos numerator, |cos|cos denominator

allocate(tot_atoms(f_ns),F3_avg(f_ns))
tot_atoms(:) = 0
F3_avg(:) = 0

do i=1,f_ns ! Iterate through species of interest (probably just OW)
    do j=1,n_f_ws(i) ! Iterate through atoms of that species
        if (pos(cart,list_f_ws(i,j)).ge.f_zmin.and.pos(cart,list_f_ws(i,j)).le.f_zmax) then
            ! Count how many atoms of interest are in the Z-range
            tot_atoms(i) = tot_atoms(i) + 1
            ! If atom is in Z-region of interest, calculate it's first coordination shell
            first_coord_shell(:,:) = 0
            size_first_coord_shell = 0
            do k=1,n_f_ws(i) ! Iterate through other atoms of species
                if (j/=k.and.pos(cart,list_f_ws(i,k)).ge.(f_zmin-f_cut).and.pos(cart,list_f_ws(i,k)).le.(f_zmax+f_cut) then
                    dx = pos(1,list_f_ws(i,k))-pos(1,list_f_ws(i,j))
                    dy = pos(2,list_f_ws(i,k))-pos(2,list_f_ws(i,j))
                    dz = pos(3,list_f_ws(i,k))-pos(3,list_f_ws(i,j))
                    call images(cart,0,1,1,icell,dx,dy,dz)
                    dsq = dx*dx + dy*dy + dz*dz
                    if (dsq <= f_cut*f_cut) then
                        size_first_coord_shell = size_first_coord_shell + 1
                        ! if size_first_coord_shell > 10, or whatever allocated size, warn & stop
                        if (size_first_coord_shell > 10) then
                            write(99,*) "WARNING: (F3) first coordination shell for atom ", list_f_ws(i,j), &
                                        ", at frame ", count, " exceeds 10 atoms!"
                            EXIT
                        first_coord_shell(size_first_coord_shell,1) = list_f_ws(i,k)
                        first_coord_shell(size_first_coord_shell,2) = dx
                        first_coord_shell(size_first_coord_shell,3) = dy
                        first_coord_shell(size_first_coord_shell,4) = dz
                        first_coord_shell(size_first_coord_shell,5) = dsq
                    endif
                endif
            enddo
            F3_atom = 0
            do k=1,size_first_coord_shell-1
                do l=k+1,size_first_coord_shell
                    ! Calculate F3 for k-j-l
                    kdotl = first_coord_shell(k,2)*first_coord_shell(l,2) + &
                            first_coord_shell(k,3)*first_coord_shell(l,3) + &
                            first_coord_shell(k,4)*first_coord_shell(l,4)
                    cos2_num = kdotl*abs(kdotl)
                    cos2_den = first_coord_shell(k,5)*first_coord_shell(l,5)
                    F3_part = (cos2_num/cos2_den + 0.1111)**2
                    ! Add F3/#combinations to total F3
                    F3_atom = F3_atom + 2*F3_part/(size_first_coord_shell**2 - size_first_coord_shell)
                enddo
            enddo
            ! f3 with a double loop on size_first_coord_shell
            ! this is where you color
            ! Calculate average F3 for the frame (per species)
            F3_avg(i) = F3_avg(i) + F3_atom
            ! later on - clustering
        endif
    enddo
    F3_avg(i) = F3_avg(i)/tot_atoms
enddo

! color by f3 value
            

end subroutine clathrates_f3

end module MOD_rings
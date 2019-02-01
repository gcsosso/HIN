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

end module MOD_rings
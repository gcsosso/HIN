module MOD_bondorder

contains

! Creates output directories and files for F3 and/or F4 calculations 
subroutine bondorder_alloc(ql)

    implicit none

    integer :: ql
    
    write(99,*) "We are calculating the q(I2) order parameter." ql
    open(unit=236, file='hin_structure.out.q(I2)_order', status='unknown')

end subroutine bondorder_alloc

subroutine bondorder(ql,q_zmin,q_zmax,q_cut,n_f_ow,list_f_ow,counter &
                      time,cart,icell,pos,nat,natformat)

    implicit none

    real :: q_zmin, q_zmax, q_cut
    real :: time
    integer :: i, j, cart, nat
    integer :: counter                      ! Frame
    integer :: n_f_ow                       ! Number of OW atoms
    real :: icell(cart*cart)
    integer, allocatable :: list_f_ow(:)    ! Atom indices of OW
    integer :: tot_atoms                    ! Count of number of atoms for which F3 is calculated
    real :: first_coord_shell(20,4)      ! First coordination shell of the current atom: (dx, dy, dz, dsq)
    integer :: first_coord_shell_ndx(20)    ! First coordination shell atom indices
    integer :: size_first_coord_shell       ! Size of first coordination shell
    real, allocatable :: pos(:,:)
    character*100 :: natformat
    real :: qlb_mol(n_f_ow), w_oz(n_f_ow)
    real :: qlb_atom

    tot_atoms = 0
    w_oz(:) = 0.0

    do i=1,n_f_ow ! Iterate through OW atoms
        if (pos(cart,list_f_ow(i)).ge.q_zmin.and.pos(cart,list_f_ow(i)).le.q_zmax) then
            ! Count how many atoms of interest are in the Z-range
            tot_atoms = tot_atoms + 1
            w_oz(tot_atoms) = pos(cart,list_f_ow(i))
            
            ! If atom is in Z-region of interest, calculate it's first coordination shell
            call compute_clath_coord_shell(i,first_coord_shell,first_coord_shell_ndx,size_first_coord_shell,n_f_ow, &
                                           f_zmin,f_zmax,q_cut,cart,icell,counter,pos,list_f_ow)
            
            ! Compute the F4 parameter for the atom
            call compute_qlb(i,ql,qlb_atom,first_coord_shell,first_coord_shell_ndx,size_first_coord_shell)
                
            qlb_mol(tot_atoms) = qlb_atom
        endif
    enddo

    write(n_mol_format,*) 2*tot_atoms
    
    write(236,'('//adjustl(n_mol_format)//'F11.4)') (w_oz(i), i=1,tot_atoms)
    
    write(236,'('//adjustl(n_mol_format)//'F11.4)') (qlb_mol(i), i=1,tot_atoms)

end subroutine bondorder


! Computes first coordination shell
subroutine compute_clath_coord_shell(i,first_coord_shell,first_coord_shell_ndx,size_first_coord_shell,n_f_ow, &
                                     q_zmin,q_zmax,q_cut,cart,icell,counter,pos,list_f_ow)

    implicit none
    
    integer :: i, j, cart, counter          ! Atom numbers for central OW, other OW atoms
    integer :: n_f_ow                       ! Number of OW atoms
    integer, allocatable :: list_f_ow(:)    ! Atom indices of OW
    integer, allocatable :: tot_atoms(:)    ! Count of number of atoms for which F3 is calculated
    real :: first_coord_shell(20,4)         ! First coordination shell of the current atom: (dx, dy, dz, dsq)
    integer :: first_coord_shell_ndx(20)    ! First coordination shell atom indices
    integer :: size_first_coord_shell       ! Size of first coordination shell
    real :: dx, dy, dz                      ! X, Y and Z distances between two atoms
    real :: dsq                             ! Square distance between two atoms
    real :: icell(cart*cart)
    real :: q_zmin, q_zmax, q_cut
    real, allocatable :: pos(:,:)
    
    size_first_coord_shell = 0
    do j=1,n_f_ow ! Iterate through other atoms of species
        if (i/=j.and.pos(cart,list_f_ow(j)).ge.(q_zmin-q_cut).and.pos(cart,list_f_ow(j)).le.(q_zmax+q_cut)) then
            dx = pos(1,list_f_ow(j))-pos(1,list_f_ow(i))
            dy = pos(2,list_f_ow(j))-pos(2,list_f_ow(i))
            dz = pos(3,list_f_ow(j))-pos(3,list_f_ow(i))
            call images(cart,0,1,1,icell,dx,dy,dz)
            dsq = dx*dx + dy*dy + dz*dz
            if (dsq <= q_cut*q_cut) then
                size_first_coord_shell = size_first_coord_shell + 1
                ! if size_first_coord_shell > 20, or whatever allocated size, warn & stop
                if (size_first_coord_shell >= 20) then
                    write(99,*) "WARNING: (F3) first coordination shell for atom ", list_f_ow(i), &
                                ", at frame ", counter, " exceeds 20 atoms!"
                    EXIT
                endif
                first_coord_shell_ndx(size_first_coord_shell) = list_f_ow(j)
                first_coord_shell(size_first_coord_shell,1) = dx
                first_coord_shell(size_first_coord_shell,2) = dy
                first_coord_shell(size_first_coord_shell,3) = dz
                first_coord_shell(size_first_coord_shell,4) = dsq
            endif
        endif
    enddo

end subroutine compute_clath_coord_shell


! Computes F3 order parameter for an atom
subroutine compute_qlb(i,ql,qlb_atom,first_coord_shell,first_coord_shell_ndx,size_first_coord_shell)

    implicit none
    
    integer :: i, m, ql
    real :: first_coord_shell(20,4)         ! First coordination shell of the current atom: (dx, dy, dz, dsq)
    integer :: size_first_coord_shell       ! Size of first coordination shell
    integer :: first_coord_shell_ndx(20)    ! First coordination shell atom indices
    real :: qlb_atom                        ! ql(i) bar parameter for atom i
    real :: j_dot_k, cos2_num, cos2_den     ! j.k, |cos|cos numerator, |cos|cos denominator
    real :: sigma
    complex :: qlmb
    real, parameter :: Pi = 3.14159
    
    sigma = 0.0

    if (size_first_coord_shell.gt.0) then 
    qlb_atom = 0
    do m=-ql,ql
        call compute_qlmb(i,m,qlmb)
        sigma = sigma + real(qlmb)**2+aimag(qlmb)**2
    enddo

    qlb_atom = sqrt(4*Pi*sigma/(2*ql+1))
    
    else
        qlb_atom = -2
    end if

end subroutine compute_qlb

subroutine compute_qlmb(i,ql,m,qlmb,first_coord_shell,first_coord_shell_ndx,size_first_coord_shell)
    
    implicit none

    integer :: i, j, m, ql
    real :: first_coord_shell(20,4)
    integer :: size_first_coord_shell
    integer :: first_coord_shell_ndx(20)    ! First coordination shell atom indices
    complex :: qlmb, qlm, sigma
    
    call compute_qlm(i,qlm)
    sigma = qlm

    do j=1,size_first_coord_shell
        call compute_qlm(j,qlm)
        sigma = sigma + qlm
    enddo

    qlmb = sigma/size_first_coord_shell

end subroutine compute_qlmb

subroutine compute_qlm(j,qlm)
    
    implicit none

    integer :: j
    complex :: qlm

end subroutine compute_qlm

end module MOD_bondorder

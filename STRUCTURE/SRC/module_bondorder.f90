module MOD_bondorder

contains

! Creates output directories and files for F3 and/or F4 calculations 
subroutine bondorder_alloc(l)

    implicit none

    integer :: l
	 character*100 :: fileloc
    
    write(99, '(A,I0.2,A)') 'We are calculating the q', l, ' order parameter.'
	 write(fileloc, '(A,I0.2,A)') 'hin_structure.out.q', l, '_order'
	 
    open(unit=236, file=fileloc, status='unknown')

end subroutine bondorder_alloc

subroutine bondorder(l,q_zmin,q_zmax,q_cut,counter,list_f_ow,n_f_ow, &
                      time,cart,icell,pos,nat,natformat,sym)

    implicit none

    real :: q_zmin, q_zmax, q_cut
    real :: time
    integer :: i, ii, l, cart, nat
    integer :: counter                      ! Frame
    real :: icell(cart*cart)
    integer :: tot_atoms                    ! Count of number of atoms for which F3 is calculated
    real :: first_coord_shell(20,4)      ! First coordination shell of the current atom: (dx, dy, dz, dsq)
    integer :: first_coord_shell_ndx(20)    ! First coordination shell atom indices
    integer :: size_first_coord_shell       ! Size of first coordination shell
    real, allocatable :: pos(:,:)
    character*100 :: natformat
    real :: ql_mol(n_f_ow), qlb_mol(n_f_ow), w_oz(n_f_ow)
    real :: ql_atom, qlb_atom
	 character*100 :: n_mol_format
	 character*4, allocatable :: sym(:)
	 integer, allocatable :: list_f_ow(:)
	 integer :: n_f_ow

    tot_atoms = 0
    w_oz(:) = 0.0
	 ql_mol(:) = 0.0
	 qlb_mol(:) = 0.0

    do ii=1,n_f_ow ! Iterate through atoms
	 	  i = list_f_ow(ii)
        if (pos(cart,i).ge.q_zmin.and.pos(cart,i).le.q_zmax) then
            ! Count how many atoms of interest are in the Z-range
            tot_atoms = tot_atoms + 1
            w_oz(tot_atoms) = pos(cart,i)
            
            ! If atom is in Z-region of interest, calculate it's first coordination shell
            call compute_first_coord_shell(ii,first_coord_shell,first_coord_shell_ndx,size_first_coord_shell, &
                                           q_zmin,q_zmax,q_cut,cart,icell,counter,pos,n_f_ow,list_f_ow,sym)
            
            ! Compute the local ql parameter for the atom
            call compute_ql(ii,l,ql_atom,first_coord_shell,first_coord_shell_ndx,size_first_coord_shell, &
										cart,icell,q_zmin,q_zmax,q_cut,pos,counter,n_f_ow,list_f_ow,sym)
										
            ! Compute the averaged ql parameter for the atom
            call compute_qlb(ii,l,qlb_atom,first_coord_shell,first_coord_shell_ndx,size_first_coord_shell, &
										cart,icell,q_zmin,q_zmax,q_cut,pos,counter,n_f_ow,list_f_ow,sym)
                
            ql_mol(tot_atoms) = ql_atom
            qlb_mol(tot_atoms) = qlb_atom
        endif
    enddo

    write(n_mol_format,*) tot_atoms
    
    write(236,'('//adjustl(n_mol_format)//'F11.4)') (w_oz(i), i=1,tot_atoms)
    write(236,'('//adjustl(n_mol_format)//'F11.4)') (ql_mol(i), i=1,tot_atoms)
    write(236,'('//adjustl(n_mol_format)//'F11.4)') (qlb_mol(i), i=1,tot_atoms)

end subroutine bondorder


! Computes first coordination shell
subroutine compute_first_coord_shell(ii,first_coord_shell,first_coord_shell_ndx,size_first_coord_shell, &
                                     q_zmin,q_zmax,q_cut,cart,icell,counter,pos,n_f_ow,list_f_ow,sym)

    implicit none
    
    integer :: ii, i, jj, j, cart, counter, n_f_ow          ! Atom numbers for central OW, other OW atoms
    integer, allocatable :: tot_atoms(:)    ! Count of number of atoms for which F3 is calculated
    real :: first_coord_shell(20,4)         ! First coordination shell of the current atom: (dx, dy, dz, dsq)
    integer :: first_coord_shell_ndx(20)    ! First coordination shell atom indices
    integer :: size_first_coord_shell       ! Size of first coordination shell
    real :: dx, dy, dz                      ! X, Y and Z distances between two atoms
    real :: dsq                             ! Square distance between two atoms
    real :: icell(cart*cart)
    real :: q_zmin, q_zmax, q_cut
    real, allocatable :: pos(:,:)
	 character*4, allocatable :: sym(:)
	 integer, allocatable :: list_f_ow(:)
    
    size_first_coord_shell = 0
    do jj=1,n_f_ow ! Iterate through other atoms
	     i = list_f_ow(ii)
	 	  j = list_f_ow(jj)
        if (ii/=jj.and.pos(cart,j).ge.(q_zmin-2*q_cut).and.pos(cart,j).le.(q_zmax+2*q_cut)) then
            dx = pos(1,j)-pos(1,i)
            dy = pos(2,j)-pos(2,i)
            dz = pos(3,j)-pos(3,i)
            call images(cart,0,1,1,icell,dx,dy,dz)
            dsq = dx*dx + dy*dy + dz*dz
            if (dsq <= q_cut*q_cut) then
                size_first_coord_shell = size_first_coord_shell + 1
                ! if size_first_coord_shell > 20, or whatever allocated size, warn & stop
                if (size_first_coord_shell >= 20) then
                    write(99,*) "WARNING: (F3) first coordination shell for atom ", i, &
                                ", at frame ", counter, " exceeds 20 atoms!"
                    EXIT
                endif
                first_coord_shell_ndx(size_first_coord_shell) = jj
                first_coord_shell(size_first_coord_shell,1) = dx
                first_coord_shell(size_first_coord_shell,2) = dy
                first_coord_shell(size_first_coord_shell,3) = dz
                first_coord_shell(size_first_coord_shell,4) = dsq
            endif
        endif
    enddo

end subroutine compute_first_coord_shell


! Computes local ql order parameter for an atom
subroutine compute_ql(ii,l,ql_atom,first_coord_shell,first_coord_shell_ndx,size_first_coord_shell, &
								cart,icell,q_zmin,q_zmax,q_cut,pos,counter,n_f_ow,list_f_ow,sym)

    implicit none
    
    integer :: ii, m, l, counter
    real :: first_coord_shell(20,4)         ! First coordination shell of the current atom: (dx, dy, dz, dsq)
    integer :: size_first_coord_shell       ! Size of first coordination shell
    integer :: first_coord_shell_ndx(20)    ! First coordination shell atom indices
    real :: ql_atom                         ! ql(i) parameter for atom i
    real :: j_dot_k, cos2_num, cos2_den     ! j.k, |cos|cos numerator, |cos|cos denominator
    real :: sigma
    complex :: qlm
    real, parameter :: Pi = 3.14159
	 integer :: cart, n_f_ow
    real :: icell(cart*cart)
    real :: q_zmin, q_zmax, q_cut
    real, allocatable :: pos(:,:)
	 character*4, allocatable :: sym(:)
	 integer, allocatable :: list_f_ow(:)
    
    sigma = 0.0
	 
    if (size_first_coord_shell.gt.0) then 
    	  ql_atom = 0
    	  do m=-l,l
        		call compute_qlm(ii,l,m,qlm,cart,icell,q_zmin,q_zmax,q_cut,pos,counter,n_f_ow,list_f_ow,sym)
        		sigma = sigma + real(qlm)**2+aimag(qlm)**2
    	  enddo

    	  ql_atom = sqrt(4*Pi*sigma/(2*l+1))
    
    else
        ql_atom = -2
    end if

end subroutine compute_ql


! Computes averaged ql order parameter for an atom
subroutine compute_qlb(ii,l,qlb_atom,first_coord_shell,first_coord_shell_ndx,size_first_coord_shell, &
								cart,icell,q_zmin,q_zmax,q_cut,pos,counter,n_f_ow,list_f_ow,sym)

    implicit none
    
    integer :: ii, m, l, counter
    real :: first_coord_shell(20,4)         ! First coordination shell of the current atom: (dx, dy, dz, dsq)
    integer :: size_first_coord_shell       ! Size of first coordination shell
    integer :: first_coord_shell_ndx(20)    ! First coordination shell atom indices
    real :: qlb_atom                        ! ql(i) bar parameter for atom i
    real :: j_dot_k, cos2_num, cos2_den     ! j.k, |cos|cos numerator, |cos|cos denominator
    real :: sigma
    complex :: qlmb
    real, parameter :: Pi = 3.14159
	 integer :: cart, n_f_ow
    real :: icell(cart*cart)
    real :: q_zmin, q_zmax, q_cut
    real, allocatable :: pos(:,:)
	 character*4, allocatable :: sym(:)
	 integer, allocatable :: list_f_ow(:)
    
    sigma = 0.0
	 
    if (size_first_coord_shell.gt.0) then 
    	  qlb_atom = 0
    	  do m=-l,l
        		call compute_qlmb(ii,l,m,qlmb,first_coord_shell,first_coord_shell_ndx,size_first_coord_shell, &
								cart,icell,q_zmin,q_zmax,q_cut,pos,counter,n_f_ow,list_f_ow,sym)
        		sigma = sigma + real(qlmb)**2+aimag(qlmb)**2
    	  enddo

    	  qlb_atom = sqrt(4*Pi*sigma/(2*l+1))
    
    else
        qlb_atom = -2
    end if

end subroutine compute_qlb

subroutine compute_qlmb(ii,l,m,qlmb,first_coord_shell,first_coord_shell_ndx,size_first_coord_shell, &
								cart,icell,q_zmin,q_zmax,q_cut,pos,counter,n_f_ow,list_f_ow,sym)
    
    implicit none

    integer :: ii, fi, m, l, counter
    real :: first_coord_shell(20,4)
    integer :: size_first_coord_shell
    integer :: first_coord_shell_ndx(20)    ! First coordination shell atom indices
    complex :: qlmb, qlm, sigma
	 integer :: cart, n_f_ow
    real :: icell(cart*cart)
    real :: q_zmin, q_zmax, q_cut
    real, allocatable :: pos(:,:)
	 character*4, allocatable :: sym(:)
	 integer, allocatable :: list_f_ow(:)
    
    call compute_qlm(ii,l,m,qlm,cart,icell,q_zmin,q_zmax,q_cut,pos,counter,n_f_ow,list_f_ow,sym)
    sigma = qlm

    do fi=1,size_first_coord_shell
        call compute_qlm(first_coord_shell_ndx(fi),l,m,qlm,cart,icell,q_zmin,q_zmax,q_cut,pos,counter,n_f_ow,list_f_ow,sym)
        sigma = sigma + qlm
    enddo

    qlmb = sigma/(size_first_coord_shell+1)

end subroutine compute_qlmb

subroutine compute_qlm(ii,l,m,qlm,cart,icell,q_zmin,q_zmax,q_cut,pos,counter,n_f_ow,list_f_ow,sym)
    
    implicit none

    integer :: ii, fj, l, m, counter
    complex :: qlm, Ylm, sigma
	 real :: first_coord_shell(20,4)
	 integer :: size_first_coord_shell
	 integer :: first_coord_shell_ndx(20)
	 integer :: cart, n_f_ow
    real :: icell(cart*cart)
    real :: q_zmin, q_zmax, q_cut
    real, allocatable :: pos(:,:)
	 character*4, allocatable :: sym(:)
	 integer, allocatable :: list_f_ow(:)
	 real :: th, ph
	 
	 call compute_first_coord_shell(ii,first_coord_shell,first_coord_shell_ndx,size_first_coord_shell, &
                                   q_zmin,q_zmax,q_cut,cart,icell,counter,pos,n_f_ow,list_f_ow,sym)
	 
	 sigma = (0.0, 0.0)
	 do fj=1,size_first_coord_shell
	 		if (first_coord_shell(fj,4).eq.0.0) then ; th = 0.0
			else ; th = acos(first_coord_shell(fj,3)/sqrt(first_coord_shell(fj,4))) ; end if
	 		if (first_coord_shell(fj,1).eq.0.0) then ; ph = 0.0
			else ; ph = atan(first_coord_shell(fj,2)/first_coord_shell(fj,1)) ; end if
	 		call compute_Ylm(Ylm,l,m,th,ph)
	 		sigma = sigma + Ylm
	 enddo
	 
	 if (size_first_coord_shell.ne.0) then
	 		qlm = sigma/size_first_coord_shell
	 else
	 		qlm = (0.0, 0.0)
	 end if

end subroutine compute_qlm

subroutine compute_Ylm(Ylm,l,m,th,ph)
	 
	 implicit none
	 
	 integer :: l, m
	 complex :: Ylm
	 real :: th, ph
    real, parameter :: Pi = 3.14159
	 
	 if (l.eq.3) then
	 		if (m.eq.-3) then
				Ylm = 0.125*sqrt(35.0/Pi)*(sin(th)**3)*exp(-3.0*cmplx(0,1)*ph)
			else if (m.eq.-2) then
				Ylm = 0.25*sqrt(52.5/Pi)*(sin(th)**2)*(cos(th))*exp(-2.0*cmplx(0,1)*ph)
			else if (m.eq.-1) then
				Ylm = 0.125*sqrt(21.0/Pi)*(sin(th))*(5.0*(cos(th)**2)-1.0)*exp(-cmplx(0,1)*ph)
			else if (m.eq.0) then
				Ylm = 0.25*sqrt(7.0/Pi)*(5.0*(cos(th)**3)-3.0*cos(th))
			else if (m.eq.1) then
				Ylm = -0.125*sqrt(21.0/Pi)*(sin(th))*(5.0*(cos(th)**2)-1.0)*exp(cmplx(0,1)*ph)
			else if (m.eq.2) then
				Ylm = 0.25*sqrt(52.5/Pi)*(sin(th)**2)*(cos(th))*exp(2.0*cmplx(0,1)*ph)
			else
				Ylm = -0.125*sqrt(35.0/Pi)*(sin(th)**3)*exp(3.0*cmplx(0,1)*ph)
			end if
	 else if (l.eq.4) then
	 		if (m.eq.-4) then
				Ylm = 0.1875*sqrt(17.5/Pi)*(sin(th)**4)*exp(-4.0*cmplx(0,1)*ph)
			else if (m.eq.-3) then
				Ylm = 0.375*sqrt(35.0/Pi)*(sin(th)**3)*cos(th)*exp(-3.0*cmplx(0,1)*ph)
			else if (m.eq.-2) then
				Ylm = 0.375*sqrt(2.5/Pi)*(sin(th)**2)*(7.0*(cos(th)**2)-1.0)*exp(-2.0*cmplx(0,1)*ph)
			else if (m.eq.-1) then
				Ylm = 0.375*sqrt(5.0/Pi)*(sin(th))*(7.0*(cos(th)**3)-3.0*cos(th))*exp(-cmplx(0,1)*ph)
			else if (m.eq.0) then
				Ylm = 0.1875*sqrt(1.0/Pi)*(35.0*(cos(th)**4)-30.0*(cos(th)**2)+3.0)
			else if (m.eq.1) then
				Ylm = -0.375*sqrt(5.0/Pi)*(sin(th))*(7.0*(cos(th)**3)-3.0*cos(th))*exp(cmplx(0,1)*ph)
			else if (m.eq.2) then
				Ylm = 0.375*sqrt(2.5/Pi)*(sin(th)**2)*(7.0*(cos(th)**2)-1.0)*exp(2.0*cmplx(0,1)*ph)
			else if (m.eq.3) then
				Ylm = -0.375*sqrt(35.0/Pi)*(sin(th)**3)*(cos(th))*exp(3.0*cmplx(0,1)*ph)
			else
				Ylm = 0.1875*sqrt(17.5/Pi)*(sin(th)**4)*exp(4.0*cmplx(0,1)*ph)
			end if
	 else if (l.eq.6) then
	 		if (m.eq.-6) then
				Ylm = 0.015625*sqrt(3003.0/Pi)*(sin(th)**6)*exp(-6.0*cmplx(0,1)*ph)
			else if (m.eq.-5) then
				Ylm = 0.09375*sqrt(1001.0/Pi)*(sin(th)**5)*cos(th)*exp(-5.0*cmplx(0,1)*ph)
	 		else if (m.eq.-4) then
				Ylm = 0.09375*sqrt(45.5/Pi)*(sin(th)**4)*(11.0*(cos(th)**2)-1.0)*exp(-4.0*cmplx(0,1)*ph)
			else if (m.eq.-3) then
				Ylm = 0.03125*sqrt(1365.0/Pi)*(sin(th)**3)*(11.0*(cos(th)**3)-3.0*cos(th))*exp(-3.0*cmplx(0,1)*ph)
			else if (m.eq.-2) then
				Ylm = 0.015625*sqrt(1365.0/Pi)*(sin(th)**2)*(33.0*(cos(th)**4)-18.0*(cos(th)**2)+1.0)*exp(-2.0*cmplx(0,1)*ph)
			else if (m.eq.-1) then
				Ylm = 0.0625*sqrt(136.5/Pi)*(sin(th))*(33.0*(cos(th)**5)-30.0*(cos(th)**3)+5.0*cos(th))*exp(-cmplx(0,1)*ph)
			else if (m.eq.0) then
				Ylm = 0.03125*sqrt(13.0/Pi)*(231.0*(cos(th)**6)-315.0*(cos(th)**4)+105.0*(cos(th)**2)-5.0)
			else if (m.eq.1) then
				Ylm = -0.0625*sqrt(136.5/Pi)*(sin(th))*(33.0*(cos(th)**5)-30.0*(cos(th)**3)+5.0*cos(th))*exp(cmplx(0,1)*ph)
			else if (m.eq.2) then
				Ylm = 0.015625*sqrt(1365.0/Pi)*(sin(th)**2)*(33.0*(cos(th)**4)-18.0*(cos(th)**2)+1.0)*exp(2.0*cmplx(0,1)*ph)
			else if (m.eq.3) then
				Ylm = -0.03125*sqrt(1365.0/Pi)*(sin(th)**3)*(11.0*(cos(th)**3)-3.0*cos(th))*exp(3.0*cmplx(0,1)*ph)
			else if (m.eq.4) then
				Ylm = 0.09375*sqrt(45.5/Pi)*(sin(th)**4)*(11.0*(cos(th)**2)-1.0)*exp(4.0*cmplx(0,1)*ph)
			else if (m.eq.5) then
				Ylm = -0.09375*sqrt(1001.0/Pi)*(sin(th)**5)*cos(th)*exp(5.0*cmplx(0,1)*ph)
			else
				Ylm = 0.015625*sqrt(3003.0/Pi)*(sin(th)**6)*exp(6.0*cmplx(0,1)*ph)
			end if
	 else
	 		Ylm = sqrt(((2*l+1)/(4*Pi))*(factorial(l-m)/factorial(l+m)))*compute_Plm(l,m,cos(th))*exp(cmplx(0,1)*m*ph)
	 end if
	 
end subroutine compute_Ylm

recursive function compute_Plm(l,m,x) result(Plm)

	 implicit none
	 
	 integer :: l, m
	 real :: x, Plm
	 
	 if (m.lt.0) then
	 		Plm = ((-1)**(-m))*(factorial(l+m)/factorial(l-m))*compute_Plm(l,-m,x)
	 else if (m.eq.l) then
	 		Plm = ((-1)**m)*factorial2(2*m-1)*((1-x**2)**(m/2))
	 else if ((m+1).eq.l) then
	 		Plm = x*(2*m+1)*compute_Plm(m,m,x)
	 else
	 		Plm = (x*(2*l-1)*compute_Plm(l-1,m,x) - (l+m-1)*compute_Plm(l-2,m,x))/(l-m)
	 end if
	 
end function compute_Plm

recursive function factorial(n) result(nf)

	 implicit none
	 
	 integer, intent(in) :: n
	 integer :: nf
	 
	 if (n.le.1) then
	 		nf = 1
	 else
	 		nf = n*factorial(n-1)
	 end if
	 
end function factorial

recursive function factorial2(n) result(nf)

	 implicit none
	 
	 integer, intent(in) :: n
	 integer :: nf
	 
	 if (n.le.1) then
	 		nf = 1
	 else
	 		nf = n*factorial(n-2)
	 end if
	 
end function factorial2

end module MOD_bondorder

module MOD_bondorder

contains

! Creates output directories and files for F3 and/or F4 calculations 
subroutine bondorder_alloc(l)

   implicit none

   integer :: l
   character*100 :: fileloc

   write(99, '(A,I0.1,A)') 'We are calculating the q', l, ' order parameter.'
   write(fileloc, '(A,I0.1,A)') 'hin_structure.out.q', l, '_order'

   open(unit=240+l, file=fileloc, status='unknown')
   !NB if more l values are allowed, file indices need to be rethinked

end subroutine bondorder_alloc

subroutine bondorder_t4_alloc()

   implicit none

   write(99, '(A,I0.1,A)') 'We are calculating the t4 order parameter.'

   open(unit=250, file='hin_structure.out.t4', status='unknown')
   open(unit=251, file='hin_structure.out.t4.color', status='unknown')

end subroutine bondorder_t4_alloc

subroutine bondorder(l, q_cut, qd_cut, qt_cut, counter, list_filtered, n_filtered, filt_param, switch_filt_param, max_shell, &
                     time, cart, icell, pos, nat, natformat, sym, switch_ql, switch_qd, switch_qt, switch_t4, qlb_io)

   implicit none
	 
   integer, parameter :: dp = kind(1.d0)
   real :: q_cut, qd_cut, qt_cut
   real :: time
   integer :: i, ii, l, cart, nat, m, max_shell
   integer :: counter                           ! Frame
   real :: icell(cart*cart)
   real :: first_coord_shell(max_shell,4)       ! First coordination shell of the current atom: (dx, dy, dz, dsq)
   integer :: first_coord_shell_ndx(max_shell)  ! First coordination shell atom indices
   integer :: size_first_coord_shell            ! Size of first coordination shell
   real :: pos(:,:), filt_param(:)
   character*100 :: natformat
   integer :: n_filtered(2)
   real(dp) :: ql_mol(n_filtered(1)), qlb_mol(n_filtered(1)), qlt_mol(n_filtered(1))
   real(dp) :: t4_mol(n_filtered(1)), t4_col(nat)
   integer :: size_first(n_filtered(2))
   real(dp) :: qlb_io(:)
   real(dp) :: ql_atom, qlb_atom, qlt_atom
   character*100 :: n_mol_format
   character*4 :: sym(:)
   logical(1) :: switch_ql, switch_qd, switch_qt, switch_t4, switch_filt_param
   integer :: list_filtered(:,:)
   complex(dp) :: qlm_all(-l:l,n_filtered(2))
   real, parameter :: Pi = 3.14159
   
   ql_mol(:) = 0.0
   qlb_mol(:) = 0.0
   qlt_mol(:) = 0.0
   t4_col(:) = 0.0
   
   do ii=1,n_filtered(2) ; i = list_filtered(2,ii)
      call compute_first_coord_shell(list_filtered(2,ii),first_coord_shell,first_coord_shell_ndx,size_first_coord_shell, &
                                 q_cut,cart,icell,counter,pos,n_filtered(2),list_filtered(2,:),sym,max_shell)
      
      if (size_first_coord_shell.gt.0) then
         do m=-l,l
            call compute_qlm(ii,l,m,qlm_all(m,ii),cart,icell,q_cut,pos,counter,sym, &
                                  first_coord_shell,size_first_coord_shell,max_shell)
         end do
      else
      
        qlm_all(:,ii) = (0.0, 0.0)
      end if
      size_first(ii) = size_first_coord_shell
   end do

   do ii=1,n_filtered(1) ; i = list_filtered(1,ii)

      ! Compute the ql parameter for the atom
      if (switch_ql) then
         call compute_ql(ii,l,ql_atom,size_first(ii),qlm_all,n_filtered(2))
         ql_mol(ii) = ql_atom
      end if
      
      ! Compute the averaged ql parameter for the atom
      if (switch_qd) then
         call compute_first_coord_shell(list_filtered(1,ii),first_coord_shell,first_coord_shell_ndx,size_first_coord_shell, &
                                        qd_cut,cart,icell,counter,pos,n_filtered(2),list_filtered(2,:),sym,max_shell)
         call compute_qlb(ii,l,qlb_atom,first_coord_shell_ndx,size_first_coord_shell,qlm_all, &
                          size_first(ii),max_shell,n_filtered(2))
         qlb_mol(ii) = qlb_atom
      end if

      ! Compute the Tianshu ql parameter for the atom
      if (switch_qt) then
         call compute_first_coord_shell(list_filtered(1,ii),first_coord_shell,first_coord_shell_ndx,size_first_coord_shell, &
                                        qt_cut,cart,icell,counter,pos,n_filtered(2),list_filtered(2,:),sym,max_shell)
         call compute_qlt(ii,l,qlt_atom,first_coord_shell_ndx,size_first_coord_shell,qlm_all,size_first(ii), &
                          max_shell,n_filtered(2))
         qlt_mol(ii) = qlt_atom
      end if
      
      ! Compute t4 if appropriate
      if (switch_t4.and.(l.eq.6)) then
         if (qlb_io(ii).ne.-2) then
            t4_mol(ii) = (1/Pi)*acos((qlb_io(ii)-0.36)/&
                                sqrt((qlb_io(ii)-0.36)**2+(qlb_mol(ii)-0.16)**2))
         else
            t4_mol(ii) = -2
         end if
         t4_col(i) = t4_mol(ii)
      end if
   end do

   if (switch_t4.and.(l.eq.3)) qlb_io = qlb_mol
   write(n_mol_format,*) n_filtered(1)
   write(240+l,'('//adjustl(n_mol_format)//'I7)') (list_filtered(1,ii), ii=1,n_filtered(1))
   if (switch_filt_param) write(240+l,'('//adjustl(n_mol_format)//'F11.4)') (filt_param(ii), ii=1,n_filtered(1))
   if (switch_ql) write(240+l,'('//adjustl(n_mol_format)//'F12.8)') (ql_mol(ii), ii=1,n_filtered(1))
   if (switch_qd) write(240+l,'('//adjustl(n_mol_format)//'F12.8)') (qlb_mol(ii), ii=1,n_filtered(1))
   if (switch_qt) write(240+l,'('//adjustl(n_mol_format)//'F12.8)') (qlt_mol(ii), ii=1,n_filtered(1))
   if (switch_t4.and.(l.eq.6)) then
      write(250,'('//adjustl(n_mol_format)//'I7)') (list_filtered(1,ii), ii=1,n_filtered(1))
      if (switch_filt_param) write(250,'('//adjustl(n_mol_format)//'F11.4)') (filt_param(ii), ii=1,n_filtered(1))
      write(250,'('//adjustl(n_mol_format)//'F12.8)') (t4_mol(ii), ii=1,n_filtered(1))
      write(251,'('//adjustl(natformat)//'F11.4)') (t4_col(i), i=1,nat)
   end if
   
end subroutine bondorder


! Computes first coordination shell
subroutine compute_first_coord_shell(i,first_coord_shell,first_coord_shell_ndx,size_first_coord_shell, &
                                     q_cut,cart,icell,counter,pos,n_filtered,list_filtered,sym,max_shell)

   implicit none

   integer :: i, jj, j, cart, counter, n_filtered  ! Atom numbers for central OW, other OW atoms
   real :: first_coord_shell(max_shell,4)          ! First coordination shell of the current atom: (dx, dy, dz, dsq)
   integer :: first_coord_shell_ndx(max_shell)     ! First coordination shell atom indices
   integer :: size_first_coord_shell               ! Size of first coordination shell
   real :: dx, dy, dz                              ! X, Y and Z distances between two atoms
   real :: dsq                                     ! Square distance between two atoms
   real :: icell(cart*cart)
   real :: q_cut
   real :: pos(:,:)
   character(4) :: sym(:)
   integer :: list_filtered(:)
   integer :: max_shell
   
   first_coord_shell(:,:) = 0.0
   first_coord_shell_ndx(:) = 0
   size_first_coord_shell = 0

   size_first_coord_shell = 0
   
   do jj=1,n_filtered ; j = list_filtered(jj)
      if (i.ne.j) then
         dx = pos(1,j)-pos(1,i)
         dy = pos(2,j)-pos(2,i)
         dz = pos(3,j)-pos(3,i)
         call images(cart,0,1,1,icell,dx,dy,dz)
         dsq = dx*dx + dy*dy + dz*dz
         if (dsq <= q_cut*q_cut) then
             size_first_coord_shell = size_first_coord_shell + 1
             ! if size_first_coord_shell > max_shell, or whatever allocated size, warn & stop
             if (size_first_coord_shell >= max_shell) then
                 write(99,*) "WARNING: Q first coordination shell for atom ", i, &
                             ", at frame ", counter, " exceeds ", max_shell, " atoms!"
                 EXIT
             end if
             first_coord_shell_ndx(size_first_coord_shell) = jj
             first_coord_shell(size_first_coord_shell,1) = dx
             first_coord_shell(size_first_coord_shell,2) = dy
             first_coord_shell(size_first_coord_shell,3) = dz
             first_coord_shell(size_first_coord_shell,4) = dsq
         end if
      end if
   end do

end subroutine compute_first_coord_shell


! Computes local ql order parameter for an atom
subroutine compute_ql(ii,l,ql_atom,size_first,qlm_all,n_filtered)

   implicit none

   integer, parameter :: dp = kind(1.d0)
   integer :: ii, m, l, counter, n_filtered
   integer :: size_first       				! Size of first coordination shell
   real(dp) :: ql_atom                      ! ql(i) parameter for atom i
   real :: sigma
   real, parameter :: Pi = 3.14159
   complex(dp) :: qlm_all(-l:l,n_filtered)

   sigma = 0.0

   if (size_first.gt.0) then 
      ql_atom = 0
      do m=-l,l ; sigma = sigma + real(qlm_all(m,ii))**2+aimag(qlm_all(m,ii))**2 ; end do
      ql_atom = sqrt(4.0*Pi*sigma/(2.0*l+1.0))
   else ; ql_atom = -2.0 ; end if

end subroutine compute_ql


! Computes averaged ql order parameter for an atom (Dellago)
subroutine compute_qlb(ii,l,qlb_atom,first_coord_shell_ndx,size_first_coord_shell,qlm_all,size_first,max_shell,n_filtered)

   implicit none

   integer, parameter :: dp = kind(1.d0)
   integer :: ii, m, l, max_shell, n_filtered
   integer :: size_first_coord_shell, size_first      ! Size of coordination shells
   integer :: first_coord_shell_ndx(max_shell)        ! First coordination shell atom indices
   real(dp) :: qlb_atom                               ! ql(i) bar parameter for atom i
   real :: sigma
   complex :: qlmb
   real, parameter :: Pi = 3.14159
   complex(dp) :: qlm_all(-l:l,n_filtered)

   sigma = 0.0

   if (size_first.gt.0) then 
      qlb_atom = 0
      do m=-l,l
         call compute_qlmb(ii,l,m,qlmb,first_coord_shell_ndx,size_first_coord_shell,qlm_all,max_shell,n_filtered)
         sigma = sigma + real(qlmb)**2+aimag(qlmb)**2
      end do
      qlb_atom = sqrt(4*Pi*sigma/(2*l+1))
   else ; qlb_atom = -2 ; end if

end subroutine compute_qlb


! Computes Tianshu version of ql order parameter for an atom
subroutine compute_qlt(ii,l,qlt_atom,first_coord_shell_ndx,size_first_coord_shell,qlm_all,size_first,max_shell,n_filtered)

   implicit none

   integer, parameter :: dp = kind(1.d0)
   integer :: ii, fj, jj, m, l, max_shell, n_filtered
   integer :: size_first_coord_shell, size_first   ! Size of coordination shells
   integer :: first_coord_shell_ndx(max_shell)     ! First coordination shell atom indices
   complex(dp) :: qlm_all(-l:l,n_filtered)
   complex :: qi_dot_qj, sigma
   real :: qi_sq, qj_sq
   real(dp) :: qlt_atom

   qi_sq = 0.0
   sigma = (0.0, 0.0)

   if (size_first.gt.0) then 
      do fj=1,size_first_coord_shell
         jj = first_coord_shell_ndx(fj)
         qi_dot_qj = (0.0, 0.0)
         qj_sq = 0.0
         do m=-l,l
            qi_dot_qj = qi_dot_qj + qlm_all(m,ii)*conjg(qlm_all(m,jj))
            if (fj.eq.1) qi_sq = qi_sq + real(qlm_all(m,ii))**2 + aimag(qlm_all(m,ii))**2
            qj_sq = qj_sq + real(qlm_all(m,jj))**2 + aimag(qlm_all(m,jj))**2
         end do
         sigma = sigma + qi_dot_qj/sqrt(qj_sq)
      end do
      qlt_atom = real(sigma)/(size_first_coord_shell*sqrt(qi_sq))
   else ; qlt_atom = -2.0 ; end if

end subroutine compute_qlt


subroutine compute_qlmb(ii,l,m,qlmb,first_coord_shell_ndx,size_first_coord_shell,qlm_all,max_shell,n_filtered)
    
   implicit none

   integer, parameter :: dp = kind(1.d0)
   integer :: ii, fi, m, l, max_shell, n_filtered
   integer :: size_first_coord_shell
   integer :: first_coord_shell_ndx(max_shell)    ! First coordination shell atom indices
   complex :: qlmb, qlm, sigma
   complex(dp) :: qlm_all(-l:l,n_filtered)

   qlm = qlm_all(m,ii)
   sigma = qlm

   do fi=1,size_first_coord_shell
      qlm = qlm_all(m,first_coord_shell_ndx(fi))
      sigma = sigma + qlm
   end do

   qlmb = sigma/(size_first_coord_shell+1)

end subroutine compute_qlmb

subroutine compute_qlm(ii,l,m,qlm,cart,icell,q_cut,pos,counter,sym, &
								first_coord_shell,size_first_coord_shell,max_shell)
    
   implicit none

   integer, parameter :: dp = kind(1.d0)
   integer :: ii, fj, l, m, counter, max_shell
   complex(dp) :: Ylm, sigma
   complex(dp) :: qlm
   real :: first_coord_shell(max_shell,4)
   integer :: size_first_coord_shell
   integer :: cart
   real :: icell(cart*cart)
   real :: q_cut
   real :: pos(:,:)
   character*4 :: sym(:)
   real(dp) :: th, ph
   real, parameter :: Pi = 3.14159
   
   sigma = (0.0, 0.0)
   do fj=1,size_first_coord_shell
      if (first_coord_shell(fj,4).eq.0.0) then ; th = 0.0
      else ; th = acos(first_coord_shell(fj,3)/sqrt(first_coord_shell(fj,4))) ; end if

      if (first_coord_shell(fj,1).eq.0.0) then
         if (first_coord_shell(fj,2).eq.0.0) then ; ph = 0.0
         else if (first_coord_shell(fj,2).gt.0.0) then ; ph = Pi/2
         else ; ph = 3*Pi/2 ; end if
      else if (first_coord_shell(fj,1).ge.0) then
         if (first_coord_shell(fj,2).ge.0) then ; ph = atan(first_coord_shell(fj,2)/first_coord_shell(fj,1))
         else ; ph = 3*Pi/2 - atan(first_coord_shell(fj,1)/first_coord_shell(fj,2)) ; end if
      else
         if (first_coord_shell(fj,2).gt.0) then ; ph = Pi/2 - atan(first_coord_shell(fj,1)/first_coord_shell(fj,2))
         else ; ph = Pi + atan(first_coord_shell(fj,2)/first_coord_shell(fj,1)) ; end if
      end if

      call compute_Ylm(Ylm,l,m,th,ph)
      sigma = sigma + Ylm
   end do
   
   qlm = sigma/size_first_coord_shell
   
end subroutine compute_qlm

subroutine compute_Ylm(Ylm,l,m,th,ph)
	 
   implicit none

   integer, parameter :: dp = kind(1.d0)
   integer :: l, m
   complex(dp) :: Ylm
   real(dp) :: th, ph
   real, parameter :: Pi = 3.14159

   if (l.eq.3) then
      if (m.eq.-3) then ; Ylm = 0.125*sqrt(35.0/Pi)*(sin(th)**3)*exp(-3.0*cmplx(0,1)*ph)
      else if (m.eq.-2) then ; Ylm = 0.25*sqrt(52.5/Pi)*(sin(th)**2)*(cos(th))*exp(-2.0*cmplx(0,1)*ph)
      else if (m.eq.-1) then ; Ylm = 0.125*sqrt(21.0/Pi)*(sin(th))*(5.0*(cos(th)**2)-1.0)*exp(-cmplx(0,1)*ph)
      else if (m.eq.0) then ; Ylm = 0.25*sqrt(7.0/Pi)*(5.0*(cos(th)**3)-3.0*cos(th))
      else if (m.eq.1) then ; Ylm = -0.125*sqrt(21.0/Pi)*(sin(th))*(5.0*(cos(th)**2)-1.0)*exp(cmplx(0,1)*ph)
      else if (m.eq.2) then ; Ylm = 0.25*sqrt(52.5/Pi)*(sin(th)**2)*(cos(th))*exp(2.0*cmplx(0,1)*ph)
      else ; Ylm = -0.125*sqrt(35.0/Pi)*(sin(th)**3)*exp(3.0*cmplx(0,1)*ph) ; end if
   else if (l.eq.4) then
      if (m.eq.-4) then ; Ylm = 0.1875*sqrt(17.5/Pi)*(sin(th)**4)*exp(-4.0*cmplx(0,1)*ph)
      else if (m.eq.-3) then ; Ylm = 0.375*sqrt(35.0/Pi)*(sin(th)**3)*cos(th)*exp(-3.0*cmplx(0,1)*ph)
      else if (m.eq.-2) then ; Ylm = 0.375*sqrt(2.5/Pi)*(sin(th)**2)*(7.0*(cos(th)**2)-1.0)*exp(-2.0*cmplx(0,1)*ph)
      else if (m.eq.-1) then ; Ylm = 0.375*sqrt(5.0/Pi)*(sin(th))*(7.0*(cos(th)**3)-3.0*cos(th))*exp(-cmplx(0,1)*ph)
      else if (m.eq.0) then ; Ylm = 0.1875*sqrt(1.0/Pi)*(35.0*(cos(th)**4)-30.0*(cos(th)**2)+3.0)
      else if (m.eq.1) then ; Ylm = -0.375*sqrt(5.0/Pi)*(sin(th))*(7.0*(cos(th)**3)-3.0*cos(th))*exp(cmplx(0,1)*ph)
      else if (m.eq.2) then ; Ylm = 0.375*sqrt(2.5/Pi)*(sin(th)**2)*(7.0*(cos(th)**2)-1.0)*exp(2.0*cmplx(0,1)*ph)
      else if (m.eq.3) then ; Ylm = -0.375*sqrt(35.0/Pi)*(sin(th)**3)*(cos(th))*exp(3.0*cmplx(0,1)*ph)
      else ; Ylm = 0.1875*sqrt(17.5/Pi)*(sin(th)**4)*exp(4.0*cmplx(0,1)*ph) ; end if
   else if (l.eq.6) then
      if (m.eq.-6) then ; Ylm = 0.015625*sqrt(3003.0/Pi)*(sin(th)**6)*exp(-6.0*cmplx(0,1)*ph)
      else if (m.eq.-5) then ; Ylm = 0.09375*sqrt(1001.0/Pi)*(sin(th)**5)*cos(th)*exp(-5.0*cmplx(0,1)*ph)
      else if (m.eq.-4) then ; Ylm = 0.09375*sqrt(45.5/Pi)*(sin(th)**4)*(11.0*(cos(th)**2)-1.0)*exp(-4.0*cmplx(0,1)*ph)
      else if (m.eq.-3) then
          Ylm = 0.03125*sqrt(1365.0/Pi)*(sin(th)**3)*(11.0*(cos(th)**3)-3.0*cos(th))*exp(-3.0*cmplx(0,1)*ph)
      else if (m.eq.-2) then
          Ylm = 0.015625*sqrt(1365.0/Pi)*(sin(th)**2)*(33.0*(cos(th)**4)-18.0*(cos(th)**2)+1.0)*exp(-2.0*cmplx(0,1)*ph)
      else if (m.eq.-1) then
          Ylm = 0.0625*sqrt(136.5/Pi)*(sin(th))*(33.0*(cos(th)**5)-30.0*(cos(th)**3)+5.0*cos(th))*exp(-cmplx(0,1)*ph)
      else if (m.eq.0) then ; Ylm = 0.03125*sqrt(13.0/Pi)*(231.0*(cos(th)**6)-315.0*(cos(th)**4)+105.0*(cos(th)**2)-5.0)
      else if (m.eq.1) then
          Ylm = -0.0625*sqrt(136.5/Pi)*(sin(th))*(33.0*(cos(th)**5)-30.0*(cos(th)**3)+5.0*cos(th))*exp(cmplx(0,1)*ph)
      else if (m.eq.2) then
          Ylm = 0.015625*sqrt(1365.0/Pi)*(sin(th)**2)*(33.0*(cos(th)**4)-18.0*(cos(th)**2)+1.0)*exp(2.0*cmplx(0,1)*ph)
      else if (m.eq.3) then
          Ylm = -0.03125*sqrt(1365.0/Pi)*(sin(th)**3)*(11.0*(cos(th)**3)-3.0*cos(th))*exp(3.0*cmplx(0,1)*ph)
      else if (m.eq.4) then ; Ylm = 0.09375*sqrt(45.5/Pi)*(sin(th)**4)*(11.0*(cos(th)**2)-1.0)*exp(4.0*cmplx(0,1)*ph)
      else if (m.eq.5) then ; Ylm = -0.09375*sqrt(1001.0/Pi)*(sin(th)**5)*cos(th)*exp(5.0*cmplx(0,1)*ph)
      else ; Ylm = 0.015625*sqrt(3003.0/Pi)*(sin(th)**6)*exp(6.0*cmplx(0,1)*ph) ; end if
   else
      Ylm = sqrt(((2*l+1)/(4*Pi))*(factorial(l-m)/factorial(l+m)))*compute_Plm(l,m,cos(th))*exp(cmplx(0,1)*m*ph)
   end if
	 
end subroutine compute_Ylm

recursive function compute_Plm(l,m,x) result(Plm) ! This code isn't currently used and isn't tested!

   implicit none

   integer, parameter :: dp = kind(1.d0)
   integer :: l, m
   real(dp) :: x, Plm

   if (m.lt.0) then ; Plm = ((-1)**(-m))*(factorial(l+m)/factorial(l-m))*compute_Plm(l,-m,x)
   else if (m.eq.l) then ; Plm = ((-1)**m)*factorial2(2*m-1)*((1-x**2)**(m/2))
   else if ((m+1).eq.l) then ; Plm = x*(2*m+1)*compute_Plm(m,m,x)
   else ; Plm = (x*(2*l-1)*compute_Plm(l-1,m,x) - (l+m-1)*compute_Plm(l-2,m,x))/(l-m) ; end if
	 
end function compute_Plm

recursive function factorial(n) result(nf)

   implicit none

   integer, intent(in) :: n
   integer :: nf

   if (n.le.1) then ; nf = 1
   else ; nf = n*factorial(n-1) ; end if
	 
end function factorial

recursive function factorial2(n) result(nf)

   implicit none

   integer, intent(in) :: n
   integer :: nf

   if (n.le.1) then ; nf = 1
   else ; nf = n*factorial(n-2) ; end if
	 
end function factorial2

end module MOD_bondorder

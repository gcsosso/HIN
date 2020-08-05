module MOD_clathrates

contains

! Creates output directories and files for F3 and/or F4 calculations 
subroutine clathrates_alloc(switch_f, switch_f_cls)

    implicit none

    logical(1) :: switch_f(3:4), switch_f_cls

    ! Make the tmp dir and open output files
    !call system("rm -r -f data-c ; mkdir data-c")
    open(unit=200, file='hin_structure.out.f.stats', status='unknown')
    open(unit=234, file='hin_structure.out.f_order', status='unknown')
    if (switch_f(3)) then
        write(99,*) "We are calculating the clathrate F3 order parameter."
        !open(unit=231, file='hin_structure.out.f.f3.color', status='unknown')
        !open(unit=230, file='hin_structure.out.f.f3', status='unknown')
        
        if (switch_f(4)) then
            write(99,*) "We are also calculating the clathrate F4 order parameter."
            !open(unit=241, file='hin_structure.out.f.f4.color', status='unknown')
            !open(unit=240, file='hin_structure.out.f.f4', status='unknown')
            if (switch_f_cls) then
                write(99,*) "We are also calculating clustering of ice and clathrate."
                open(unit=201, file='hin_structure.out.f.ice.patch', status='unknown')
                open(unit=202, file='hin_structure.out.f.ice.patch.color', status='unknown')
                open(unit=203, file='hin_structure.out.f.clathrate.patch', status='unknown')
                open(unit=204, file='hin_structure.out.f.clathrate.patch.color', status='unknown')
            end if
            
            write(200,*) "# Time [ps] | Average F3 | Average F4 "
        else
            write(234,*) "Z F3"
            write(200,*) "# Time [ps] | Average F3 "
        end if
    else
        write(99,*) "We are calculating the clathrate F4 order parameter."
        open(unit=241, file='hin_structure.out.f.f4.color', status='unknown')
        !open(unit=240, file='hin_structure.out.f.f4', status='unknown')
        write(234,*) "Z F4"
        write(200,*) "# Time [ps] | Average F4 "
    end if

end subroutine clathrates_alloc

subroutine clathrates(switch_f, f_cut, n_filtered, list_filtered, filt_param, switch_filt_param, counter, time, cart, &
                      icell, pos, nat, natformat, f_zbins, switch_f_cls, f3_imax, f3_cmax, f4_imax, f4_cmin)

   implicit none

   integer, parameter :: dp = kind(1.d0)
   logical(1) :: switch_f(3:4), switch_f_cls, switch_filt_param
   real :: f_cut, f3_imax, f3_cmax, f4_imax, f4_cmin
   real :: time
   integer :: i, ii, j, cart, nat
   integer :: counter                      ! Frame
   integer :: n_filtered(2)                       ! Number of OW atoms
   real :: icell(cart*cart)
   integer :: list_filtered(:,:)    ! Atom indices of OW
   real :: first_coord_shell(20,4)      ! First coordination shell of the current atom: (dx, dy, dz, dsq)
   integer :: first_coord_shell_ndx(20)    ! First coordination shell atom indices
   integer :: size_first_coord_shell       ! Size of first coordination shell
   real(dp) :: F3_atom, F4_atom                ! F3 parameter for triples, atoms
   real :: F3_avg, F4_avg                  ! F3 parameter for frame-wide avg
   real :: pos(:,:), filt_param(:)
   real :: F3_col(nat), F4_col(nat)
   character*100 :: natformat, n_mol_format
   integer :: f_zbins
   integer :: F3_zbin_len(f_zbins), F4_zbin_len(f_zbins)
   real :: F3_zbin(f_zbins,nat), F4_zbin(f_zbins,nat)
   real(dp) :: F3_mol(n_filtered(1)), F4_mol(n_filtered(1))
     
   F3_avg = 0
   F4_avg = 0
   F3_col(:) = 0.0
   F4_col(:) = 0.0
   F3_mol(:) = 0.0
   F4_mol(:) = 0.0
   F3_zbin_len(:) = 0
   F4_zbin_len(:) = 0

   do ii=1,n_filtered(1) ; i = list_filtered(1,ii)

      ! If atom is in Z-region of interest, calculate it's first coordination shell
      call compute_clath_coord_shell(i, first_coord_shell, first_coord_shell_ndx, size_first_coord_shell, n_filtered(2), &
                                     f_cut, cart, icell, counter, pos, list_filtered(2,:))
      
      if (switch_f(3)) then
          ! Compute the F3 parameter for the atom
          call compute_f3(F3_atom,first_coord_shell,size_first_coord_shell)

          ! this is where you color
          ! Calculate average F3 for the frame (per species)
          !F3_avg = F3_avg + F3_atom
          !F3_col(list_f_ow(i)) = F3_atom
          F3_mol(ii) = F3_atom
          !do j=1,f_zbins
          !    if ((pos(cart,list_f_ow(i))>=(((f_zmax-f_zmin)*(j-1)/f_zbins)+f_zmin)) .and. &
          !        (pos(cart,list_f_ow(i))<=(((f_zmax-f_zmin)*j/f_zbins)+f_zmin))) then
          !        F3_zbin_len(j) = F3_zbin_len(j) + 1
          !        F3_zbin(j,F3_zbin_len(j)) = F3_atom
          !        exit
          !    end if
          !end do
          ! later on - clustering
      end if
      if (switch_f(4)) then
          ! Compute the F4 parameter for the atom
          call compute_f4(i, F4_atom, first_coord_shell, first_coord_shell_ndx, size_first_coord_shell, &
                          cart, icell, pos)

          ! Calculate average F4 for the frame (per species)
          !F4_avg = F4_avg + F4_atom
          !F4_col(list_f_ow(i)) = F4_atom
          F4_mol(ii) = F4_atom
          !do j=1,f_zbins
          !    if (pos(cart,list_f_ow(i))<=(((f_zmax-f_zmin)*j/f_zbins)+f_zmin)) then
          !        F4_zbin_len(j) = F4_zbin_len(j) + 1
          !        F4_zbin(j,F4_zbin_len(j)) = F4_atom
          !        exit
          !    end if
          !end do
     end if
   end do
   
   !if (tot_atoms>0) then
   !    F3_avg = F3_avg/tot_atoms
   !    F4_avg = F4_avg/tot_atoms
   !end if

   write(n_mol_format,*) n_filtered(1)
   
   write(234,'('//adjustl(n_mol_format)//'I7)') (list_filtered(1,i), i=1,n_filtered(1))
   if (switch_filt_param) write(234,'('//adjustl(n_mol_format)//'F11.4)') (filt_param(i), i=1,n_filtered(1))

   if (switch_f(3)) then
     write(234,'('//adjustl(n_mol_format)//'F12.8)') (F3_mol(ii), ii=1,n_filtered(1))
   end if
   if (switch_f(4)) then
     write(234,'('//adjustl(n_mol_format)//'F12.8)') (F4_mol(ii), ii=1,n_filtered(1))
   end if
   
   !if (switch_f(3).and.switch_f(4)) then
   !  ! If we are calculating both order parameters, write to output files
   !  write(200,'(1E12.6,2(X,F12.7))') time, F3_avg, F4_avg
   !  ! Write line to color files
   !  write(231,'('//adjustl(natformat)//'F11.4)') (F3_col(i), i=1,nat)
   !  write(241,'('//adjustl(natformat)//'F11.4)') (F4_col(i), i=1,nat)
   !  ! Write lines to binned output files
   !  do j=1,f_zbins
   !      write(230,'('//adjustl(natformat)//'F11.4)') (F3_zbin(j,i), i=1,F3_zbin_len(j))
   !      write(240,'('//adjustl(natformat)//'F11.4)') (F4_zbin(j,i), i=1,F4_zbin_len(j))
   !  end do
   !else if (switch_f(3)) then
   !  ! If we are calculating only F3, write to output files
   !  write(200,'(1E12.6,X,F12.7)') time, F3_avg
   !  ! Write line to color file
   !  write(231,'('//adjustl(natformat)//'F11.4)') (F3_col(i), i=1,nat)
   !  ! Write lines to binned output files
   !  do j=1,f_zbins
   !      write(230,'('//adjustl(natformat)//'F11.4)') (F3_zbin(j,i), i=1,F3_zbin_len(j))
   !  end do
   !else
   !  ! If we are calculating only F4, write to output files
   !  write(200,'(1E12.6,X,F12.7)') time, F4_avg
   !  ! Write line to color file
   !  write(241,'('//adjustl(natformat)//'F11.4)') (F4_col(i), i=1,nat)
   !  ! Write lines to binned output files
   !  do j=1,f_zbins
   !      write(240,'('//adjustl(natformat)//'F11.4)') (F4_zbin(j,i), i=1,F4_zbin_len(j))
   !  end do
   !end if

   ! Clustering

   if (switch_f_cls) then

     if (.not.(switch_f(3).and.switch_f(4))) then
         write(99,*) "Clathrate clustering requires both F3 and F4!"
     else
         call f_clustering(F3_mol, F4_mol, nat, n_filtered(1), list_filtered(1,:), 0.0, f3_imax, -1.0, f4_imax, f_cut, &
                           pos, icell, 6, 201, 202, natformat, time) ! Cluster icy molecules
         call f_clustering(F3_mol, F4_mol, nat, n_filtered(1), list_filtered(1,:), 0.0, f3_cmax, f4_cmin, 1.0, f_cut, &
                           pos, icell, 5, 203, 204, natformat, time) ! Cluster clathrate-like molecules
     end if

   end if

end subroutine clathrates


! Computes first coordination shell
subroutine compute_clath_coord_shell(i, first_coord_shell, first_coord_shell_ndx, size_first_coord_shell, n_filtered, &
                                     f_cut, cart, icell, counter, pos, list_filtered)

   implicit none

   integer :: i, j, jj, cart, counter          ! Atom numbers for central OW, other OW atoms
   integer :: n_filtered                   ! Number of OW atoms
   integer :: list_filtered(:)             ! Atom indices of OW
   real :: first_coord_shell(20,4)         ! First coordination shell of the current atom: (dx, dy, dz, dsq)
   integer :: first_coord_shell_ndx(20)    ! First coordination shell atom indices
   integer :: size_first_coord_shell       ! Size of first coordination shell
   real :: dx, dy, dz                      ! X, Y and Z distances between two atoms
   real :: dsq                             ! Square distance between two atoms
   real :: icell(cart*cart)
   real :: f_cut
   real :: pos(:,:)

   size_first_coord_shell = 0
   do jj=1,n_filtered ; j = list_filtered(jj) ; if (i.ne.j) then
      dx = pos(1,j)-pos(1,i)
      dy = pos(2,j)-pos(2,i)
      dz = pos(3,j)-pos(3,i)
      call images(cart,0,1,1,icell,dx,dy,dz)
      dsq = dx*dx + dy*dy + dz*dz
      if (dsq <= f_cut*f_cut) then
         size_first_coord_shell = size_first_coord_shell + 1
         ! if size_first_coord_shell > 20, or whatever allocated size, warn & stop
         if (size_first_coord_shell >= 20) then
            write(99,*) "WARNING: (F3) first coordination shell for atom ", i, &
                        ", at frame ", counter, " exceeds 20 atoms!"
            EXIT
         end if
         first_coord_shell_ndx(size_first_coord_shell) = j
         first_coord_shell(size_first_coord_shell,1) = dx
         first_coord_shell(size_first_coord_shell,2) = dy
         first_coord_shell(size_first_coord_shell,3) = dz
         first_coord_shell(size_first_coord_shell,4) = dsq
      end if
   end if ; end do

end subroutine compute_clath_coord_shell


! Computes F3 order parameter for an atom
subroutine compute_f3(F3_atom, first_coord_shell, size_first_coord_shell)

    implicit none
    
	 integer, parameter :: dp = kind(1.d0)
    integer :: j, k
    real :: first_coord_shell(20,4)         ! First coordination shell of the current atom: (dx, dy, dz, dsq)
    integer :: size_first_coord_shell       ! Size of first coordination shell
    real(dp) :: F3_part, F3_atom                ! F3 parameter for triples, atoms
    real(dp) :: j_dot_k, cos2_num, cos2_den     ! j.k, |cos|cos numerator, |cos|cos denominator
    
    if (size_first_coord_shell.gt.1) then 
    F3_atom = 0
    do j=1,size_first_coord_shell-1
        do k=j+1,size_first_coord_shell
            j_dot_k = first_coord_shell(j,1)*first_coord_shell(k,1) + &
            first_coord_shell(j,2)*first_coord_shell(k,2) + &
            first_coord_shell(j,3)*first_coord_shell(k,3)
            cos2_num = j_dot_k*abs(j_dot_k)
            cos2_den = first_coord_shell(j,4)*first_coord_shell(k,4)
            F3_part = (cos2_num/cos2_den + 0.1111)**2
            ! Add F3/#combinations to total F3
            F3_atom = F3_atom + F3_part
        end do
    end do
    
!    F3_atom = 2*F3_atom/(size_first_coord_shell**2 - size_first_coord_shell)
    
    else
        F3_atom = -2
    end if

end subroutine compute_f3


! Computes F4 order parameter for an atom
subroutine compute_f4(i, F4_atom, first_coord_shell, first_coord_shell_ndx, size_first_coord_shell, &
                      cart, icell, pos)

   implicit none

   integer, parameter :: dp = kind(1.d0)
   integer :: i, j, jj, cart
   real :: first_coord_shell(20,4)         ! First coordination shell of the current atom: (dx, dy, dz, dsq)
   integer :: first_coord_shell_ndx(20)    ! First coordination shell indices
   integer :: size_first_coord_shell       ! Size of first coordination shell
   real(dp) :: F4_part, F4_atom                ! F3 parameter for triples, atoms
   real :: pos(:,:)
   real :: icell(cart*cart)
   real :: dx1, dy1, dz1, dx2, dy2, dz2    ! Tmp holders for component distance vectors
   real :: h1x, h1y, h1z, h2x, h2y, h2z    ! Positional vector components for H1 and H2 (H1-O1..O2-H2, O1 is origin)
   real :: ux, uy, uz, vx, vy, vz          ! Vectors in (H1-O1..O2), (O1..O2,H2) planes, perp to intersection
   real :: lambda, mu                      ! u = (1,lambda), v=(1,mu) under bases: {h1, o2} and {h2, o2}, resp.
   real :: h1_dot_o2, h2_dot_o2, u_dot_v
   real :: cos_phi                         ! phi is the torsional angle (H1-O1..O2-H2)

   if (size_first_coord_shell.gt.0) then
      F4_atom = 0
      do jj=1,size_first_coord_shell ; j = first_coord_shell_ndx(jj)
         ! Choose Hydrogen from O1. I.e. furthest from O2.
         dx1 = pos(1,i+1)-pos(1,j)
         dy1 = pos(2,i+1)-pos(2,j)
         dz1 = pos(3,i+1)-pos(3,j)
         call images(cart,0,1,1,icell,dx1,dy1,dz1)
         dx2 = pos(1,i+2)-pos(1,j)
         dy2 = pos(2,i+2)-pos(2,j)
         dz2 = pos(3,i+2)-pos(3,j)
         call images(cart,0,1,1,icell,dx2,dy2,dz2)
         if (dx1*dx1+dy1*dy1+dz1*dz1>dx2*dx2+dy2*dy2+dz2*dz2) then
            ! Calculate h1-o1
            h1x = pos(1,i+1)-pos(1,i)
            h1y = pos(2,i+1)-pos(2,i)
            h1z = pos(3,i+1)-pos(3,i)
            call images(cart,0,1,1,icell,h1x,h1y,h1z)
         else
            ! Calculate h1-o1
            h1x = pos(1,i+2)-pos(1,i)
            h1y = pos(2,i+2)-pos(2,i)
            h1z = pos(3,i+2)-pos(3,i)
            call images(cart,0,1,1,icell,h1x,h1y,h1z)
         end if

         ! Choose Hydrogen from O2. I.e. furthest from O1.
         dx1 = pos(1,j+1)-pos(1,i)
         dy1 = pos(2,j+1)-pos(2,i)
         dz1 = pos(3,j+1)-pos(3,i)
         call images(cart,0,1,1,icell,dx1,dy1,dz1)
         dx2 = pos(1,j+2)-pos(1,i)
         dy2 = pos(2,j+2)-pos(2,i)
         dz2 = pos(3,j+2)-pos(3,i)
         call images(cart,0,1,1,icell,dx2,dy2,dz2)
         if (dx1*dx1+dy1*dy1+dz1*dz1>dx2*dx2+dy2*dy2+dz2*dz2) then
            h2x = dx1
            h2y = dy1
            h2z = dz1
         else
            h2x = dx2
            h2y = dy2
            h2z = dz2
         end if

         ! Calculate lambda, mu
         h1_dot_o2 = h1x*first_coord_shell(jj,1) + h1y*first_coord_shell(jj,2) + h1z*first_coord_shell(jj,3)
         h2_dot_o2 = h2x*first_coord_shell(jj,1) + h2y*first_coord_shell(jj,2) + h2z*first_coord_shell(jj,3)
         lambda = -h1_dot_o2/first_coord_shell(jj,4)
         mu = -h2_dot_o2/first_coord_shell(jj,4)

         ! Calculate u, v
         ux = h1x + lambda*first_coord_shell(jj,1)
         uy = h1y + lambda*first_coord_shell(jj,2)
         uz = h1z + lambda*first_coord_shell(jj,3)
         vx = h2x + mu*first_coord_shell(jj,1)
         vy = h2y + mu*first_coord_shell(jj,2)
         vz = h2z + mu*first_coord_shell(jj,3)

         ! Calculate F4 contribution
         u_dot_v = ux*vx + uy*vy + uz*vz
         cos_phi = u_dot_v/sqrt((ux*ux+uy*uy+uz*uz)*(vx*vx+vy*vy+vz*vz))
         F4_part = 4*cos_phi**3 - 3*cos_phi

         ! Add F4/#combinations to total F4
         F4_atom = F4_atom + F4_part/size_first_coord_shell
      end do
      
   else ; F4_atom = -2 ; end if

end subroutine compute_f4


subroutine f_clustering(F3_mol, F4_mol, nat, n_filtered, list_filtered, f3_min, f3_max, f4_min, f4_max, f_cut, &
                        pos, icell, crit, patch_file, col_file, natformat, time)
    
    use dfs
    
    implicit none
    
    integer, parameter :: dp = kind(1.d0)
    integer :: nat, n_filtered, i, ii, j
    integer :: list_filtered(:)
    real(dp) :: F3_mol(n_filtered), F4_mol(n_filtered)
    real :: f3_min, f3_max, f4_min, f4_max, f_cut
    real :: pos(:,:)
    real :: posi(3), posj(3), icell(9), xdf, ydf, zdf, dist
    logical(1) :: cknn
    integer :: patch_file, col_file
    character*100 :: natformat
    real :: time
    
    ! DFS stuff
    integer :: ncr, iat, jat, nnf, voltot, crit, ncrit, vol_count, patch
    integer, allocatable :: dfs_color(:), volume_crit(:)
    
    allocate(cr_list(n_filtered), lwho(n_filtered,n_filtered), dfs_color(nat), volume_crit(n_filtered))
    ncr = 0
    cr_list(:) = 0
    
    do ii=1,n_filtered ; i = list_filtered(ii)
        if (F3_mol(ii).ge.f3_min.and.F3_mol(ii).le.f3_max&
            &.and.F4_mol(ii).ge.f4_min.and.F4_mol(ii).le.f4_max) then
            ncr = ncr + 1
            cr_list(ncr) = i
        end if
    end do
    
    allocate(graph_solid_connect(ncr,20),followgraph(ncr),volume(ncr))
    allocate(neigh(ncr),predecessor(ncr))
    ! Fill the adjacency list
    neigh(:)=0
    lwho(:,:)=0
    graph_solid_connect(:,:)=0
    predecessor(:)=0
    
    do i=1,ncr
        do j=1,ncr
            if (i.ne.j) then
                iat = cr_list(i)
                jat = cr_list(j)
                posi(:)=pos(:,iat)
                posj(:)=pos(:,jat)
                call nn (posj,posi,icell,f_cut,cknn,xdf,ydf,zdf,dist)
                if (cknn) then 
                    neigh(i) = neigh(i)+1
                    graph_solid_connect(i,neigh(i)) = j
                end if
            end if
        end do
    end do
    
    count_cls=0
    volume(:)=0
    vol_count=0
    followgraph(:)=0
    
    do i=1,ncr
        black = 0
        if (followgraph(i).eq.0) then
            count_cls = count_cls + 1
            followgraph(i) = f_explore(i)
            volume(count_cls) = black
        end if
    end do
    
    nnf=count_cls
    
    ! number of clusters (not yet filtered by crit) = nnf  
    ! volume of each cluster = volume(j) up to count
    ! list of atoms indexes that constitute the clusters = lwho
    
    count_cls=0
    dfs_color(:)=0
    voltot=0
    volume_crit(:)=0
    
    do i=1,nnf
        if (volume(i).gt.crit) then
            do j=1,volume(i)
                dfs_color(lwho(i,j)) = i
            end do
            count_cls = count_cls + 1
            voltot = voltot + volume(i)
            volume_crit(count_cls) = volume(i)
        end if
    end do

    ncrit = count_cls

    ! number of clusters greater than crit = ncrit

    ! sort the volumes: the biggest is the surface itself
    call f_sort2(volume_crit,count_cls)

    ! Here is the number of atoms within the largest hexagonal patch
    patch = volume_crit(1)
    ! Write down patch statistics: time, n. of atoms in the biggest patch, n. of atoms which meet the parameter
    write(patch_file,"(1E10.4,2i10)") time, patch, ncr

    ! Write down the colors for VMD
    write(col_file,"("//adjustl(natformat)//"i10)") (dfs_color(i), i=1,nat)

    deallocate(graph_solid_connect,followgraph,volume,neigh,predecessor,lwho,cr_list,dfs_color,volume_crit)

end subroutine f_clustering

recursive function f_explore(index) result(fat)

    use dfs
    
    implicit none
    
    integer, intent(in) :: index
    integer :: j, iat
    integer :: fat


    followgraph(index) = 1
    do j=1,neigh(index)
        iat = graph_solid_connect(index, j)
        if (followgraph(iat).eq.0) then
            predecessor(iat) = index
            followgraph(iat) = f_explore(iat)
        end if
    end do
    black = black + 1
    fat=2

    ! who's who?
    lwho(count_cls,black)=cr_list(index)

end function f_explore

subroutine f_sort2(dati, n) ! Insertion sort

    integer :: n
    integer, dimension(n) :: dati
    integer :: i, j, tmp, min, pos

    do i=1,n-1
        min = dati(i)
        pos = i
        do j=i+1,n
            if (dati(j)>min) then
                min = dati(j) 
                pos = j
            end if
        end do
        tmp = dati(i) 
        dati(i) = dati(pos)
        dati(pos) = tmp
    end do

end subroutine f_sort2

end module MOD_clathrates

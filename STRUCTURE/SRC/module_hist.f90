module MOD_hist

contains

! Creates output directories and files and reads 
subroutine hist_alloc(nat, hist_centre, resname, n_hist_cs, list_hist_cs, hist_x, hist_nbins, hist_counts, switch_hist_noh, sym)

   implicit none
	
   logical(1) :: tmp_hs_ast, switch_hist_noh
   integer :: nat, n_hist_cs, hist_nbins
   integer, allocatable :: list_hist_cs(:), hist_counts(:,:)
	character*5, allocatable :: resname(:)
   character(20) :: hist_centre
	character(5) :: hist_x
   character*4, allocatable :: sym(:)

   write(99,*) 'We are calculating the histogram.'
   
	allocate(list_hist_cs(nat))
   open(unit=111, file='hin_structure.out.hist', status='unknown')
	if (hist_x.eq.'shell') then
		call read_hshell_centre(nat, hist_centre, resname, n_hist_cs, list_hist_cs, switch_hist_noh, sym)
	else if (hist_x.eq.'z') then
		write(99,*) "Slab based histogram computing is not yet implemented."
		stop
	else
		write(99,*) "hist -x value not recognised."
		stop
	end if
	
	allocate(hist_counts(2,hist_nbins))
	hist_counts(:,:) = 0

end subroutine hist_alloc

subroutine hist(nat, hist_x, n_hist_cs, list_hist_cs, resname, pos, hist_min, hist_max, hist_nbins, &
	 				 n_filtered, list_filtered, hist_counts, cart, icell)

   implicit none
	
   integer :: nat, n_hist_cs, n_filtered(2), i, j, ii, jj, cart, bin, hist_nbins
   integer, allocatable :: list_hist_cs(:), list_filtered(:,:), hist_counts(:,:)
   character*5, allocatable :: resname(:), list_all_ws(:)
	character(5) :: hist_x
	real, allocatable :: pos(:,:)
	real :: hist_min, hist_max, dx, dy, dz, dsq, tmp_dsq
	real :: icell(cart*cart)
	
	if (hist_x.eq.'shell') then
		 do i=1,n_filtered(1)
			 ii = list_filtered(1,i)
          dsq = hist_max*hist_max+1
			 do j=1,n_hist_cs
				jj = list_hist_cs(j)
				dx = pos(1,ii)-pos(1,jj)
				dy = pos(2,ii)-pos(2,jj)
				dz = pos(3,ii)-pos(3,jj)
				call images(cart,0,1,1,icell,dx,dy,dz)
				tmp_dsq = dx*dx + dy*dy + dz*dz
            if (tmp_dsq.lt.dsq) dsq = tmp_dsq
          end do
          if ((dsq.lt.(hist_max*hist_max)).and.(dsq.ge.(hist_min*hist_min))) then
            bin = ceiling((sqrt(dsq)*real(hist_nbins))/(hist_max-hist_min))
            hist_counts(1,bin) = hist_counts(1,bin) + 1
          end if
		 end do
		 do i=1,n_filtered(2)
			 ii = list_filtered(2,i)
          dsq = hist_max*hist_max+1
			 do j=1,n_hist_cs
				jj = list_hist_cs(j)
				dx = pos(1,ii)-pos(1,jj)
				dy = pos(2,ii)-pos(2,jj)
				dz = pos(3,ii)-pos(3,jj)
				call images(cart,0,1,1,icell,dx,dy,dz)
				tmp_dsq = dx*dx + dy*dy + dz*dz
            if (tmp_dsq.lt.dsq) dsq = tmp_dsq
          end do
          if ((dsq.lt.(hist_max*hist_max)).and.(dsq.ge.(hist_min*hist_min))) then
				bin = ceiling((sqrt(dsq)*real(hist_nbins))/(hist_max-hist_min))
            if (bin.eq.0) bin = 1
				hist_counts(2,bin) = hist_counts(2,bin) + 1
          end if
		 end do
	else if (hist_x.eq.'z') then
		! Slab based histogram to go here
		write(99,*) "Slab based histogram computing is not yet implemented."
		stop
	end if
   
end subroutine hist

subroutine hist_output(hist_x, hist_min, hist_max, hist_nbins, hist_counts)
   
   implicit none
   
   integer :: hist_nbins, i
   integer, allocatable :: hist_counts(:,:)
	character(5) :: hist_x
	real :: hist_min, hist_max, binstart, binend, frac, binsize
	 
   if (hist_x.eq.'shell') then
	 	write(99,*) "We have calculated the histogram from the index file as a function of distance from the specified centre."
   else if (hist_x.eq.'z') then
      write(99,*) "We have calculated the histogram from the index file as a function of z."
   end if
   
   write(111,*) "  Start   |    End    |  No. icy  | Total No. |  Fraction"
   
   binend = hist_min
   binsize = (hist_max-hist_min)/real(hist_nbins)
   do i=1,hist_nbins
      binstart = binend
      binend = binstart + binsize
      if (hist_counts(2,i).gt.0) then
        frac = real(hist_counts(1,i))/real(hist_counts(2,i))
      else
        frac = 0
      end if
      write(111,'(F11.7,X,F11.7,X,E11.4,X,E11.4,X,F11.9)') binstart, binend, real(hist_counts(1,i)), real(hist_counts(2,i)), frac
   end do

end subroutine hist_output



subroutine read_hshell_centre(nat, hist_centre, resname, n_hist_cs, list_hist_cs, switch_hist_noh, sym)

   implicit none

   logical(1) :: centre_range, switch_hist_noh
   integer :: nat, i, n_hist_cs, delim_index, centre_start, centre_end
   integer, allocatable :: list_hist_cs(:)
   character(1) :: delim=':'
   character*5, allocatable :: resname(:)
   character(20) :: hist_centre
   character*4, allocatable :: sym(:)

   ! -centre input can be provided as atom index range or residue name
   if (verify(delim,hist_centre).eq.0) then ! If colon detected then interpret as a index range
     centre_range = .true.
     delim_index = scan(hist_centre,delim)
     read(hist_centre(1:delim_index-1),*) centre_start
     read(hist_centre(delim_index+1:),*) centre_end
   else ; centre_range = .false. ; end if ! Otherwise interpret as resname

   do i=1,nat
      if (centre_range.and.(i.ge.centre_start).and.(i.le.centre_end)) then ;
        if (.not.switch_hist_noh.or.(sym(i)(1:1).ne.'H')) then
            n_hist_cs = n_hist_cs + 1
            list_hist_cs(n_hist_cs) = i
        end if
      else if ((.not.centre_range).and.(trim(adjustl(resname(i))).eq.trim(adjustl(hist_centre)))) then ;
        if (.not.switch_hist_noh.or.(sym(i)(1:1).ne.'H')) then
            n_hist_cs = n_hist_cs + 1
            list_hist_cs(n_hist_cs) = i
        end if
      end if
	end do

end subroutine read_hshell_centre

end module MOD_hist
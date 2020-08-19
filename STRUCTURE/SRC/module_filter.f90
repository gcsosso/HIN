module MOD_filter

contains

subroutine initial_filter(nat, ns, ws, n_ws, list_ws, sym, n_all_ws, list_all_ws, centre, resname, n_cs, list_cs)

   implicit none

   logical(1) :: tmp_ws_ast, centre_range
   integer :: nat, ns, i, j, tmp_ws_len, n_all_ws, n_cs, delim_index, centre_start, centre_end
   integer, allocatable :: n_ws(:), list_ws(:,:), list_all_ws(:), list_cs(:)
   character(1) :: delim=':'
   character(5) :: tmp_ws
   character(4), allocatable :: sym(:), ws(:)
   character*5, allocatable :: resname(:)
   character(7) :: centre

   allocate(list_ws(ns,nat), n_ws(ns), list_all_ws(nat), list_cs(nat))
   n_ws(:) = 0
   n_all_ws = 0
   n_cs = 0

   ! -centre input can be provided as atom index range or residue name
   if (verify(delim,centre).eq.0) then ! If colon detected then interpret as a index range
     centre_range = .true.
     delim_index = scan(centre,delim)
     read(centre(1:delim_index-1),*) centre_start
     read(centre(delim_index+1:),*) centre_end
   else ; centre_range = .false. ; end if ! Otherwise interpret as resname

   do i=1,nat
      ! Filter by ws. An * can be used at the end of a species name to select anything with the same starting characters.
      ! E.g. HW* for both HW1 and HW2.
      if (centre_range.and.(i.ge.centre_start).and.(i.le.centre_end)) then ;
        n_cs = n_cs + 1
        list_cs(n_cs) = i
      else if ((.not.centre_range).and.(trim(adjustl(resname(i))).eq.trim(adjustl(centre)))) then ;
        n_cs = n_cs + 1
        list_cs(n_cs) = i
      end if
      do j=1,ns
         tmp_ws_ast = .false.
         tmp_ws = ws(j)
         tmp_ws_len = len(trim(adjustl(tmp_ws)))
         if (tmp_ws(tmp_ws_len:tmp_ws_len).eq.'*') then ; tmp_ws_ast = .true. ; tmp_ws_len = tmp_ws_len-1 ; end if
         if (((.not.tmp_ws_ast).and.(trim(adjustl(sym(i))).eq.trim(adjustl(tmp_ws)))) &
             .or.(tmp_ws_ast.and.(sym(i)(1:tmp_ws_len).eq.tmp_ws(1:tmp_ws_len)))) then
            n_ws(j) = n_ws(j) + 1
            n_all_ws = n_all_ws + 1
            list_ws(j, n_ws(j)) = i
            list_all_ws(n_all_ws) = i
         end if
      end do
   end do

end subroutine initial_filter


subroutine frame_filter(filter, filt_min, filt_max, op_max_cut, n_all_ws, list_all_ws, n_filtered, list_filtered, sym, ns, &
                        pos, filt_param, qlb_io, n_cs, list_cs, cart, icell)

   implicit none

   integer, parameter :: dp = kind(1.d0)
   integer :: ns, i, ii, j, n_all_ws, n_filtered(2), n_cs, cart
   integer, allocatable :: list_all_ws(:), list_filtered(:,:), list_cs(:)
   real :: filt_min, filt_max, op_max_cut, icell(cart*cart)
   real, allocatable :: pos(:,:), filt_param(:)
   character(7) :: filter
   character(4), allocatable :: sym(:)
   real(dp), allocatable :: qlb_io(:)

   allocate(list_filtered(2,n_all_ws), filt_param(n_all_ws))
   n_filtered = 0

   do ii=1,n_all_ws
      i = list_all_ws(ii)
      if (trim(adjustl(filter)).eq.'none') then
         n_filtered(:) = n_all_ws
         list_filtered(1,:) = list_all_ws ; list_filtered(2,:) = list_all_ws
      else if (trim(adjustl(filter)).eq.'z') then
        call filter_z(i, filt_min, filt_max, op_max_cut, pos(3,list_all_ws(ii)), n_filtered, list_filtered, filt_param)
      else if (trim(adjustl(filter)).eq.'shell') then
        call filter_shell(i, filt_min, filt_max, op_max_cut, n_filtered, list_filtered, filt_param, n_cs, list_cs, pos, cart, icell)
      end if
   end do
   allocate(qlb_io(n_filtered(1)))

end subroutine frame_filter


subroutine filter_z(i, zmin, zmax, op_max_cut, zpos, n_filtered, list_filtered, filt_param)

   implicit none

   integer :: i, n_filtered(2)
   real :: zmin, zmax, zpos, op_max_cut
   integer :: list_filtered(:,:)
   real :: filt_param(:)

   if ((zpos.ge.zmin).and.(zpos.le.zmax)) then
      n_filtered(1) = n_filtered(1) + 1 ; n_filtered(2) = n_filtered(2) + 1
      list_filtered(1,n_filtered) = i ; list_filtered(2,n_filtered) = i
      filt_param(n_filtered) = zpos
   else if ((zpos.ge.(zmin-op_max_cut)).and.(zpos.le.(zmax+op_max_cut))) then
      n_filtered(2) = n_filtered(2) + 1
      list_filtered(2,n_filtered) = i
   end if

end subroutine filter_z


subroutine filter_shell(i, rmin, rmax, op_max_cut, n_filtered, list_filtered, filt_param, n_cs, list_cs, pos, cart, icell)

  implicit none

  ! Arguments
  integer :: i, n_filtered(2), n_cs, cart
  real :: rmin, rmax, op_max_cut, icell(cart*cart)
  real, allocatable :: pos(:,:), filt_param(:)
  integer :: list_filtered(:,:), list_cs(:)

  ! Local
  integer :: j, jj
  real :: min_dist, i_pos(cart), j_pos(cart), xdf, ydf, zdf, r_ij

  min_dist=rmax+op_max_cut
  i_pos(1)=pos(1,i) ; i_pos(2)=pos(2,i) ; i_pos(3)=pos(3,i)
  do jj=1,n_cs
    j = list_cs(jj)
    j_pos(1)=pos(1,j) ; j_pos(2)=pos(2,j) ; j_pos(3)=pos(3,j)
    xdf=i_pos(1)-j_pos(1) ; ydf=i_pos(2)-j_pos(2) ; zdf=i_pos(3)-j_pos(3)
    call images(cart,0,1,1,icell,xdf,ydf,zdf)
    r_ij=sqrt(xdf**2.0d0+ydf**2.0d0+zdf**2.0d0)
    if (r_ij.lt.min_dist) then
      min_dist = r_ij
    end if
  end do

  if ((min_dist.ge.rmin).and.(min_dist.le.rmax)) then
    n_filtered(1) = n_filtered(1) + 1 ; n_filtered(2) = n_filtered(2) + 1
    list_filtered(1,n_filtered(1)) = i ; list_filtered(2,n_filtered(2)) = i
    filt_param(n_filtered(1)) = min_dist
  else if ((min_dist.ge.(rmin-op_max_cut)).and.(min_dist.lt.(rmax+op_max_cut))) then
    n_filtered(2) = n_filtered(2) + 1
    list_filtered(2,n_filtered(2)) = i
  end if

end subroutine filter_shell

end module MOD_filter

module MOD_filter

contains

subroutine initial_filter(nat, ns, ws, n_ws, list_ws, sym, n_all_ws, list_all_ws)
   
   implicit none
   
   logical(1) :: tmp_ws_ast
   integer :: nat, ns, i, j, tmp_ws_len, n_all_ws
   integer, allocatable :: n_ws(:), list_ws(:,:), list_all_ws(:)
   character(5) :: tmp_ws
   character(4), allocatable :: sym(:), ws(:)
   
   allocate(list_ws(ns,nat), n_ws(ns), list_all_ws(nat))
   n_ws(:) = 0
   n_all_ws = 0
   
   do i=1,nat
      ! Filter by ws. An * can be used at the end of a species name to select anything with the same starting characters.
      ! E.g. HW* for both HW1 and HW2.   
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
                        pos, filt_param, qlb_io)
   
   implicit none
   
   integer, parameter :: dp = kind(1.d0)
   integer :: ns, i, ii, j, n_all_ws, n_filtered(2)
   integer, allocatable :: list_all_ws(:), list_filtered(:,:)
   real :: filt_min, filt_max, op_max_cut
   real, allocatable :: pos(:,:), filt_param(:)
   character(7) :: filter
   character(4), allocatable :: sym(:)
   real(dp), allocatable :: qlb_io(:)
   
   allocate(list_filtered(2,n_all_ws), filt_param(n_all_ws))
   n_filtered = 0
   
   do ii=1,n_all_ws
      i = list_all_ws(ii)
      if (trim(adjustl(filter)).eq.'z') then
         call filter_z(i, filt_min, filt_max, op_max_cut, pos(3,list_all_ws(ii)), n_filtered, list_filtered, filt_param)
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

end module MOD_filter
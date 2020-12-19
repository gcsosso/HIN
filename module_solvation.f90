module MOD_solvation

contains

subroutine solvation(list_s_ws,pos,ns,n_ws,current_color,icell,cart,s_rcut,current_coord,counter,ws)
implicit none
!Arguments
integer :: ns, cart, counter
real ::  icell(cart*cart)
real, allocatable :: pos(:,:)
integer, allocatable :: n_ws(:), list_s_ws(:,:)
integer, allocatable :: current_color(:), current_coord(:,:)
character*4, allocatable :: sym(:),ws(:)
! Local variables
integer :: i, j, m, col_flg
real :: xd, yd, zd,s_rcut
real :: rij_M

current_coord(:,:)=0

do i=1,ns-1 !loops over all species given in input minus OW species

    do j=1,n_ws(1) !loops over all OWs n_ws(1)= OW which is why OW
                  ! needs to be the first argument after -species instruction
        if (current_color(list_s_ws(1,j)).eq.1) then
           col_flg = 1 !icy
        else
            if (current_color(list_s_ws(1,j)).eq.0) then
            col_flg = 2 !liquid
            end if    
        end if 
        
        do m = 1,n_ws(i+1) ! if more than one of species i, loop over them
             xd = pos(1,list_s_ws(1,j)) - pos(1,list_s_ws(i+1,m))
             yd = pos(2,list_s_ws(1,j)) - pos(2,list_s_ws(i+1,m))
             zd = pos(3,list_s_ws(1,j)) - pos(3,list_s_ws(i+1,m))
  
             call images(cart,0,1,1,icell,xd,yd,zd) 
             rij_M=sqrt(xd**2.0 + yd**2.0 + zd**2.0) 

             if (rij_M.lt.s_rcut) then
                  current_coord(i,col_flg)=current_coord(i,col_flg)+1
            
            ! write(*,*)'ts, species, current_coord ,icy ,liquid , water idx,',counter, ws(i+1), current_coord(i,col_flg),current_coord(i,1),current_coord(i,2), j,xd,yd,zd
             end if
        end do 
    end do
    write(1314,*),'ts  ,species  ,icy  ,liquid', counter, ws(i+1), current_coord(i,1), current_coord(i,2)
end do

end subroutine solvation

end module MOD_solvation

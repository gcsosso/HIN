subroutine nn (posj,posi,icell,rcut,cknn,xdf,ydf,zdf,rsqdf)

implicit none

integer, parameter :: cart=3
real, intent(in) :: posj(cart), posi(cart), icell(cart*cart), rcut
logical, intent(inout) :: cknn
real :: xdf, ydf, zdf, rsqdf !(=>dist elsewhere...)

xdf=posj(1)-posi(1)
ydf=posj(2)-posi(2)
zdf=posj(3)-posi(3)
call images(cart,0,1,1,icell,xdf,ydf,zdf)
rsqdf=sqrt(xdf**2.0+ydf**2.0+zdf**2.0)
if (rsqdf.lt.rcut) then ! mk is a neighbor of l1
   cknn=.true.
else
   cknn=.false.
endif
   
return

end subroutine nn

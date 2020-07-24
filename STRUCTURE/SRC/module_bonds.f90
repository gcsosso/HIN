module MOD_bonds

contains

subroutine bonds_alloc(nz_bAVE,b_zmax,b_zmin,b_dz,pdbon_AVE,npairs,b_bins,cn_AVE,npairs_cn)

implicit none

! Local

! Arguments
integer :: b_bins, nz_bAVE, npairs, npairs_cn
real :: b_zmax, b_zmin, b_dz
real, allocatable :: pdbon_AVE(:,:,:), cn_AVE(:,:)

nz_bAVE=int((b_zmax-b_zmin)/b_dz)
if (nz_bAVE.lt.1) then
   write(99,*) "Your mesh along z (bonds section) is ill-defined. Modify b_zmin, b_zmax and/or b_dz"
   stop
endif
allocate(pdbon_AVE(nz_bAVE,npairs,b_bins),cn_AVE(nz_bAVE,npairs_cn))
pdbon_AVE(:,:,:)=0.0
cn_AVE(:,:)=0.0


end subroutine bonds_alloc

subroutine bonds(b_zmin,b_zmax,b_dz,pos,icell,pdbon,cart,ns,n_ws,list_ws, &
                 ws,b_rcut,npairs,b_bins,b_bmin,b_bmax,cn,npairs_cn,cn_AVE,pdbon_AVE)

implicit none

! Arguments
integer :: cart, ns, npairs, b_bins, npairs_cn
integer, allocatable :: n_ws(:), list_ws(:,:)
real :: b_zmin,b_zmax,b_dz, icell(cart*cart), b_bmin, b_bmax
real, allocatable :: pos(:,:), b_rcut(:), pdbon_AVE(:,:,:),cn_AVE(:,:)
character*4, allocatable :: ws(:), ck_ws(:,:)

! Local
integer :: nz, i, j, k, counter, flag, ipair, l, m, ibin, pij
integer, allocatable :: nbl(:,:)
real :: posl(cart), posm(cart), lb, ub, xdf, ydf, zdf, dist, half_dr, rstep
real :: uv(cart)
real, allocatable :: pdbon(:,:,:), norm(:,:), cn(:,:)
logical(1) :: cknn

nz=int((b_zmax-b_zmin)/b_dz)
allocate(pdbon(nz,npairs,b_bins),norm(nz,npairs),cn(nz,npairs_cn))
pdbon(:,:,:)=0.0
cn(:,:)=0.0

half_dr=(abs(b_bmax-b_bmin))/(real(b_bins)*2.0)

allocate(ck_ws(ns*ns,2))
counter=0
ipair=0
allocate(nbl(nz,ns))
nbl(:,:)=0

! Count the number of atoms for each species in every slice
do i=1,ns
   do l=1,n_ws(i)
      do k=1,nz
         lb=b_zmin+(k-1)*b_dz
         ub=lb+b_dz
         if (pos(cart,list_ws(i,l)).gt.lb.and.pos(cart,list_ws(i,l)).lt.ub) then
            nbl(k,i)=nbl(k,i)+1
         endif
      enddo
   enddo
enddo  

do i=1,ns
   do j=1,ns
      counter=counter+1
      ck_ws(counter,1)=ws(i) ; ck_ws(counter,2)=ws(j)
      flag=0
      do k=1,counter
         if (ck_ws(counter,1).eq.ck_ws(k,2).and.ck_ws(counter,2).eq.ck_ws(k,1).and.ck_ws(counter,1).ne.ck_ws(counter,2)) flag=1
      enddo
      if (flag.eq.0) then ! We have a unique pair
         ipair=ipair+1
         do l=1,n_ws(i) ! Loop over species ws(i) - l
            do m=l,n_ws(j) ! Loop over species ws(j) - m | Avoid double
                           !counting of the bonds, but cn = cn*2
               do k=1,nz ! Get the bonds for each pair along z...
                  lb=b_zmin+(k-1)*b_dz
                  ub=lb+b_dz
                  if (list_ws(i,l).eq.list_ws(j,m)) cycle
                  if (pos(cart,list_ws(i,l)).gt.lb.and.pos(cart,list_ws(i,l)).lt.ub) then 
                     ! Atom l of species i is within the z slice of interest
                     if (pos(cart,list_ws(j,m)).gt.lb.and.pos(cart,list_ws(j,m)).lt.ub) then
                        ! Atom l of species i is within the z slice of interest
                        ! Now get the bond length for these two:
                        posl(:)=pos(:,list_ws(i,l)) ; posm(:)=pos(:,list_ws(j,m))
                        call nn (posl,posm,icell,b_rcut(ipair),cknn,xdf,ydf,zdf,dist)
                        if (cknn) then
                           cn(k,ipair)=cn(k,ipair)+1.0
                        !! DEBUG
                        !! Get the ho-oh-xy_ANGLE
                        !write(987,*) asin(-zdf/dist)*57.2958
                        !! END DEBUG
                        endif
                        do ibin=1,b_bins
                           rstep=((abs(b_bmax-b_bmin)/real(b_bins))*real(ibin))+b_bmin
                           if (dist.gt.rstep-half_dr.and.dist.lt.rstep+half_dr) then
                              pdbon(k,ipair,ibin)=pdbon(k,ipair,ibin)+1.0
                           endif
                        enddo
                     endif
                  endif
               enddo                
            enddo
         enddo
      endif
   enddo
enddo

! Get the norm
counter=0
ipair=0
norm(:,:)=0.0
do i=1,ns
   do j=1,ns
      counter=counter+1
      ck_ws(counter,1)=ws(i) ; ck_ws(counter,2)=ws(j)
      flag=0
      do k=1,counter
         if (ck_ws(counter,1).eq.ck_ws(k,2).and.ck_ws(counter,2).eq.ck_ws(k,1).and.ck_ws(counter,1).ne.ck_ws(counter,2)) flag=1
      enddo
      if (flag.eq.0) then ! We have a unique pair
         ipair=ipair+1
         do k=1,nz
            do ibin=1,b_bins
               if (pdbon(k,ipair,ibin).gt.0.0) then
                  norm(k,ipair)=norm(k,ipair)+((abs(b_bmax-b_bmin)/real(b_bins))*pdbon(k,ipair,ibin)) 
               endif
            enddo
         enddo
      endif
   enddo
enddo

! Normalize pdens
counter=0
ipair=0
do i=1,ns
   do j=1,ns
      counter=counter+1
      ck_ws(counter,1)=ws(i) ; ck_ws(counter,2)=ws(j)
      flag=0
      do k=1,counter
         if (ck_ws(counter,1).eq.ck_ws(k,2).and.ck_ws(counter,2).eq.ck_ws(k,1).and.ck_ws(counter,1).ne.ck_ws(counter,2)) flag=1
      enddo
      if (flag.eq.0) then ! We have a unique pair
         ipair=ipair+1
         do k=1,nz
            do ibin=1,b_bins
               if (pdbon(k,ipair,ibin).gt.0.0) then
                  pdbon(k,ipair,ibin)=pdbon(k,ipair,ibin)/norm(k,ipair)
               else
                  pdbon(k,ipair,ibin)=0.0
               endif
            enddo
         enddo
      endif
   enddo
enddo

! Normalize coord
ipair=0
counter=0
cn(:,:)=cn(:,:)*2.0 ! No double counting in the m,l loop above...

! Normalize up to npairs
do i=1,ns
   do j=1,ns
      counter=counter+1
      ck_ws(counter,1)=ws(i) ; ck_ws(counter,2)=ws(j)
      flag=0
      do k=1,counter
         if (ck_ws(counter,1).eq.ck_ws(k,2).and.ck_ws(counter,2).eq.ck_ws(k,1).and.ck_ws(counter,1).ne.ck_ws(counter,2)) flag=1
      enddo
      if (flag.eq.0) then ! We have a unique pair
         ipair=ipair+1
         do k=1,nz
            if (nbl(k,i).gt.0) then
               cn(k,ipair)=cn(k,ipair)/real(nbl(k,i))
            else
               cn(k,ipair)=0.0
            endif
         enddo
      endif
   enddo
enddo

   
counter=0
ipair=0
pij=0!! Normalize the npairs+cn

do i=1,ns
   do j=1,ns
      counter=counter+1
      ck_ws(counter,1)=ws(i) ; ck_ws(counter,2)=ws(j)
      flag=0
      do k=1,counter
         if (ck_ws(counter,1).eq.ck_ws(k,2).and.ck_ws(counter,2).eq.ck_ws(k,1).and.ck_ws(counter,1).ne.ck_ws(counter,2)) flag=1
      enddo
      if (flag.eq.0) then ! We have a unique pair
         ipair=ipair+1
         if (ws(i).ne.ws(j)) then
            pij=pij+1
            cn(:,npairs+pij)=cn(:,ipair)
            do k=1,nz 
               if (nbl(k,i).gt.0.and.nbl(k,j).gt.0) then
                  cn(k,npairs+pij)=(cn(k,npairs+pij)*real(nbl(k,i)))/real(nbl(k,j))
               else
                  cn(k,npairs+pij)=0.0
               endif
            enddo
         endif
      endif
   enddo
enddo 

pdbon_AVE(:,:,:)=pdbon_AVE(:,:,:)+pdbon(:,:,:)
cn_AVE(:,:)=cn_AVE(:,:)+cn(:,:)

deallocate(pdbon,norm,ck_ws,nbl,cn)

return

end subroutine bonds

end module MOD_bonds

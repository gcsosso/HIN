! pair correlation functions
 implicit none
 integer i,itim,iatm,jatm,ntim,natm,nblk,idum,deltat,tstart,tend,nr,ir,&
&        ityp,jtyp,ntyp,ndim,icount
 real(8) dr,half_dr,xi,yi,zi,dx,dy,dz,rr,fact0,fact,pi4,vol,cell(3),rho,r2
 character(len=70) nfile,label!,cdum
 character(len=2) lb1,lb2
 logical tp_exist
 integer,parameter :: mxtyp=15, two=2
 integer,allocatable :: nattyp(:)
 real(8),allocatable :: gr(:,:),rx(:,:),ry(:,:),rz(:,:),rad(:),grint(:,:),tr(:),tmpsum(:), smgr(:,:)
 character(len=2),allocatable :: lb_spc(:,:),atyp(:)

 tp_exist=.false.
 do i=1,3
  call getarg(1+i,label)
  read(label,*,end=40) cell(i)
 enddo
 tp_exist=.true.
 40 continue
 if(.not. tp_exist) then
   write(6,'(a)') 'Usage:'
   call getarg(0,label)
   write(6,'(a,a)') trim(label),' <file>, A, B, C, #radial pt, deltat, tstart, [tend]'
   stop
 endif
 
 tend=0
 call getarg(5,label); read(label,*) nr
 call getarg(6,label); read(label,*) deltat
 call getarg(7,label); read(label,*) tstart
 call getarg(8,label); read(label,*,err=100,end=100) tend
 100 continue
!DEBUG
!write(6,*) (cell(i),i=1,3)
!write(6,*) deltat,tstart,tend; stop"OK"
!ENDDEBUG

!Build radial (regular) mesh, initialize g(r)
 allocate(rad(nr))
 dr=(minval(cell(:))-5.d-1)/dble(nr)
 dr=(minval(cell(:)))/dble(2*nr)
 half_dr=dr/2.d0
 do ir=1,nr
   rad(ir)=(dble(ir-1)+0.5d0)*dr
 enddo

 call getarg(1,nfile)
 open(5,file=trim(nfile),action='read')

 read(5,*) natm; nblk=natm+2; ntim=1
 20 continue
 read(5,*,end=30) label
 ntim=ntim+1
 goto 20
 30 continue
 !write(*,*) "DEBUG, ntim, nblk, mod", ntim, nblk, mod(ntim,nblk)
 if(mod(ntim,nblk)/=0) stop "mod(ntim,nblk)/=0"
 ntim=ntim/nblk-1
 !write(6,'(a,5i7)') 'ntim, nr, deltat, tstart, tend=',ntim,nr,deltat,tstart,tend !DEBUG
 rewind(5)

!Read xyz file
 allocate(rx(natm,0:ntim),ry(natm,0:ntim),rz(natm,0:ntim),lb_spc(natm,0:ntim))
 do itim=0,ntim
   !write(*,*) "reading", itim
   read(5,*) idum
   read(5,'(a)') label
!DEBUG
!  read(label,'(a,a,i8,a)') lb1,lb2,idum,cdum
!  write(6,*) idum,itim
!ENDDEBUG
   do iatm=1,natm
     read(5,*) lb_spc(iatm,itim),rx(iatm,itim),ry(iatm,itim),rz(iatm,itim)
   enddo
 enddo
 close(5)

!How many species? What?
 allocate(nattyp(0:mxtyp),atyp(mxtyp)); nattyp(:)=0; nattyp(0)=natm
!At least one atom of type lb_spc(1,0)
 atyp(1)=lb_spc(1,0); nattyp(1)=1; ntyp=1
 do iatm=2,natm

   tp_exist=.false.
   tploop: do ityp=1,ntyp
     if(lb_spc(iatm,0)==atyp(ityp))then
       nattyp(ityp)=nattyp(ityp)+1
       tp_exist=.true.
       exit tploop
     endif
   enddo tploop
   if(.not. tp_exist)then
     ntyp=ntyp+1; if(ntyp>mxtyp) stop "increase mxtyp"
     atyp(ntyp)=lb_spc(iatm,0)
     nattyp(ntyp)=1
   endif

 enddo
 ndim=ntyp**2 !ndim=((ntyp+1)*ntyp)/2
 allocate(gr(nr,0:ndim),smgr(nr,0:ndim),tr(nr),tmpsum(0:ndim),grint(nr,0:ndim)); gr(:,:)=0.d0; smgr(:,:)=0.d0
!DEBUG
 !write(6,'(a,5i7)') 'ntyp, ndim, natm=              ',ntyp,ndim,natm
 !write(6,'(10i5)') (nattyp(ityp),ityp=0,ntyp)
 !write(6,'(2x,a,10(3x,a))') 'sum',(atyp(ityp),ityp=1,ntyp)
!stop"OK"
!ENDDEBUG

 if(tend>ntim.or.tend==0) tend=ntim
 icount=0
 do itim=tstart,tend
   !write(*,*) "conf", itim
   if(mod(itim-tstart,deltat)==0)then
   icount=icount+1
   !write(*,*) itim, icount

     do iatm=1,natm
       lb1=lb_spc(iatm,itim); xi=rx(iatm,itim); yi=ry(iatm,itim); zi=rz(iatm,itim)
       do jatm=1,natm
         lb2=lb_spc(jatm,itim)
         dx=rx(jatm,itim)-xi; dx=dx-dnint(dx/cell(1))*cell(1)
         dy=ry(jatm,itim)-yi; dy=dy-dnint(dy/cell(2))*cell(2)
         dz=rz(jatm,itim)-zi; dz=dz-dnint(dz/cell(3))*cell(3)
         rr=dsqrt(dx*dx + dy*dy + dz*dz)
radloop: do ir=1,nr
!write(6,'(i4,2e16.8)') ir,rad(ir),rr
           if(rr>rad(ir)-half_dr.and.rr<=rad(ir)+half_dr)then
             gr(ir,0)=gr(ir,0)+1.d0
             do ityp=1,ntyp

               idum=ntyp*(ityp-1)
               do jtyp=1,ntyp
                 if( (lb1==atyp(jtyp).and.lb2==atyp(ityp)) )then
                   i=idum+jtyp; if(i>ndim) stop "out of gr"
                   gr(ir,i)=gr(ir,i)+1.d0
                 endif
               enddo

!!! !              Omo
!!!                if(lb1==atyp(ityp).and.lb2==atyp(ityp))then
!!!                  gr(ir,ityp)=gr(ir,ityp)+1.d0
!!!                endif
!!! !              Etero
!!!                do jtyp=ityp+1,ntyp
!!! !write(6,'(2i4,4(" "a))') ir,idum,lb1,atyp(ityp),lb2,atyp(jtyp)
!!!                  if( (lb1==atyp(ityp).and.lb2==atyp(jtyp)) .or. &
!!! &                    (lb1==atyp(jtyp).and.lb2==atyp(ityp)) )then
!!!                    idum=ntyp*ityp-ityp*(ityp+1)/2+jtyp; if(idum>ndim) stop "out of gr"
!!!                    gr(ir,idum)=gr(ir,idum)+1.d0
!!!                  endif
!!!                enddo

             enddo
             exit radloop
           endif
         enddo radloop
       enddo
     enddo

   endif
 enddo

!DEBUG
 do ir=1,nr
!write(*,*) "smooth", ir
! SMOOTH
do i=0,ndim
   smgr(ir,i)=(gr(ir-2,i)+2.0d0*gr(ir-1,i)+3.0d0*gr(ir,i)+2.0d0*gr(ir+1,i)+gr(ir+2,i))/9.0d0
   if (ir.le.two) then
      smgr(ir,i)=gr(3,i)
   endif
   if (ir.ge.nr-2) then
      smgr(ir,i)=gr(nr,i)
   endif
enddo
! 
 write(99,'(30es16.8)') rad(ir),(smgr(ir,i),i=0,ndim),(gr(ir,i),i=0,ndim)
 enddo
!ENDDEBUG

!Output files
 !total.dat : r | TOT g(r) | TOT intg(r)
 !          : PAIRS: GeGe, GeTe, ecc... same format: r | PAIR g(r) | PAIR intg(r)
 !          THESE ARE THE FILES TO BE PLOTTED !! NOW SMOOTHED !!
 open(50,file='total.dat',action='readwrite')
 ! all.dat and integr.dat are not that useful (confusion...)
 open(40,file='all.dat',action='readwrite')
 open(30,file='integr.dat',action='readwrite')
 do ityp=1,ntyp
   idum=ntyp*(ityp-1)
   do jtyp=1,ntyp
     i=idum+jtyp
     nfile=trim(atyp(jtyp))//trim(atyp(ityp))//'.dat'
     open(50+i,file=trim(nfile),action='readwrite')
   enddo
 enddo
!Normalize and output g(r)
 pi4=16.d0*datan(1.d0)
 vol=cell(1)*cell(2)*cell(3)
 rho=dble(natm)/vol
 fact0=1.d0/(pi4*rho*dble(icount)*dr)
 grint(:,:)=0.d0; tmpsum(:)=0.d0
 do ir=1,nr
   !write(*,*) "output", ir
   r2=rad(ir)**2
   grint(ir,:)=tmpsum(:)
   ! DEBUG
   gr(ir,0)=gr(ir,0)*fact0/dble(natm)/r2
   smgr(ir,0)=smgr(ir,0)*fact0/dble(natm)/r2
!a gr(ir,0)=0.d0
   do ityp=1,ntyp
     idum=ntyp*(ityp-1)
     fact=1.d0/(pi4*dble(icount)*dr*dble(nattyp(ityp))/vol)
     do jtyp=1,ntyp
       i=idum+jtyp
!!!!   gr(ir,i)=gr(ir,i)/( pi4*dble(icount)*dr *dble(nattyp(ityp)**2)*r2 /vol)
!!!!   grint(ir,i)=grint(ir,i)+pi4* (dble(nattyp(ityp))/vol) *dr*r2*gr(ir,i)
       gr(ir,i)=gr(ir,i)*fact/dble(nattyp(jtyp))/r2
       smgr(ir,i)=smgr(ir,i)*fact/dble(nattyp(jtyp))/r2
       grint(ir,i)=grint(ir,i) + pi4*r2*(dble(nattyp(jtyp))/vol)*dr* gr(ir,i)
!a     gr(ir,0)=gr(ir,0)+ gr(ir,i)*dble(nattyp(ityp))*dble(nattyp(jtyp))/dble(natm**2)
     enddo
   enddo
   grint(ir,0)=grint(ir,0) + pi4*r2*(dble(natm)/vol)*dr* gr(ir,0)
   tmpsum(:)=grint(ir,:)
   tr(ir)=pi4*dble(natm)*rad(ir)*gr(ir,0)/vol
   do i=0,ndim
     ! SM TRICK
     if (ir.lt.3) then
     smgr(ir,i)=smgr(3,i)
     endif
     write(50+i,'(30es16.8)') rad(ir),gr(ir,i),smgr(ir,i),grint(ir,i)
   enddo
   write(40,'(30es16.8)') rad(ir),gr(ir,0),tr(ir),(gr(ir,i),i=1,ndim)
 enddo

 deallocate(tr)
 allocate(tr(ntyp))

 do ir=1,nr
   tr(:)=0.d0
   do ityp=1,ntyp
     idum=ntyp*(ityp-1)
     do jtyp=1,ntyp
       i=idum+jtyp
       tr(ityp)=tr(ityp)+grint(ir,i)
     enddo
   enddo
   write(30,'(30es16.8)') rad(ir),(tr(ityp),ityp=1,ntyp)
 enddo

 do i=0,ndim
   close(50+i)
 enddo
 close(40)
 close(30)

 deallocate(rx,ry,rz,lb_spc,rad,gr,smgr,nattyp,atyp,grint,tr)
 stop
 end


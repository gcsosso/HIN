module MOD_clusters

contains

subroutine clusters_alloc(outxtc,ns,ws,hw_ex,ohstride,n_ws,list_ws,sfile,vmd_exe,n_cls_AVE)

   implicit none

   ! Local
   integer :: i, j, k, flag, tmplist
   character*1 :: ch
   character*100 :: command, command1, command2, fcommand, pstring, pstring_C

   ! Arguments
   integer :: ns, ohstride
   integer, allocatable :: n_ws(:), list_ws(:,:)
   real :: n_cls_AVE
   character*3 :: outxtc, hw_ex
   character*4, allocatable :: ws(:)
   character*100 :: sfile, vmd_exe

   ! Check whether we are writing down the .xtc(s). If not, no way, plumed need
   ! them on a frame-2-frame basis
   if (trim(adjustl(outxtc)).ne.'yes') then
      write(99,*) "PLUMED2 needs the .xtc to be written (set OUTXTC yes in hin_structure.in)"
      stop
   endif
   ! Check whether we have both oxygens and hydrogens, and that the HW_EXC is yes
   flag=0
   do i=1,ns
      if(ws(i).eq.'OW') flag=flag+1
      if(ws(i).eq.'HW1') flag=flag+1
   enddo
   if (flag.ne.2) then
      write(99,*) "In order to have clusters, you need to read both oxygens and hydrogens..."
      stop
   endif
   if (trim(adjustl(hw_ex)).ne.'yes') then
      write(99,*) "Set HW_EXC to yes, please..."
      stop
   endif
   ! Write the rigth atomic indexes in plumed.dat
   command="cp plumed.dat_TEMPLATE plumed.dat"
   call system(command)

   if (ohstride.eq.0) then
      write(99,*) "Doing indexes by hand is not implemented - yet..."
      stop
   endif
   ! oxygens
   do i=1,ns
      if (ws(i).eq.'OW') then
         do j=1,n_ws(i)
            if (j.eq.1) then
               tmplist=list_ws(i,j)
            else
               if (list_ws(i,j).ne.tmplist+ohstride) then
                  write(99,*) "Sort oxygens and hydrogens properly, mate..."
                  stop
               else
                  tmplist=list_ws(i,j)
               endif
            endif
         enddo
         write(pstring,*) list_ws(i,1), "-", list_ws(i,n_ws(i)), ":", ohstride
         pstring_C=" "
         do k=1,len_trim(pstring)
           ch=pstring(k:k)
           if (ch.ne." ") then
             pstring_C=trim(pstring_C) // ch
           endif
         end do
         pstring_C=trim(adjustl(pstring_C))
         command1="sed -i 's/TBC_OX/"
         command2="/' plumed.dat"
         fcommand=trim(command1)//trim(pstring_C)//trim(command2)
         call system(fcommand)
      endif
   enddo
   ! hydrogens
   do i=1,ns
      if (ws(i).eq.'HW1') then ! this gets HW1 AND HW2
         do j=1,n_ws(i)
            if (j.eq.1) then
               tmplist=list_ws(i,j)
            else
               if (list_ws(i,j).eq.tmplist+1) then !HW2 1 Assumes it follows straight after HW1
               else if (list_ws(i,j).eq.tmplist+ohstride) then !HW1
                  tmplist=list_ws(i,j)
               else
                  write(99,*) "Sort oxygens and hydrogens properly, mate..."
                  stop
               endif
            endif
         enddo

         write(pstring,*) list_ws(i,1), "-", list_ws(i,n_ws(i))-1, ":", ohstride, ",", &
                          list_ws(i,1)+1, "-", list_ws(i,n_ws(i)), ":", ohstride
         pstring_C=" "
         do k=1,len_trim(pstring)
           ch=pstring(k:k)
           if (ch.ne." ") then
             pstring_C=trim(pstring_C) // ch
           endif
         end do
         pstring_C=trim(adjustl(pstring_C))
         command1="sed -i 's/TBC_HY/"
         command2="/' plumed.dat"
         fcommand=trim(command1)//trim(pstring_C)//trim(command2)
         call system(fcommand)
      endif
   enddo
   ! Convert conf.gro in conf.pdb
   open(unit=69, file='tmp.tcl', status='unknown')
   write(69,*) "package require topotools"
   write(69,*) "package require psfgen"
   write(69,*) "package require pbctools"
   write(69,*) "mol load gro  ", trim(sfile)
   write(69,*) "set sel [atomselect top all]"
   write(69,*) "$sel writepdb conf.pdb"
   write(69,*) "quit"
   close(69)
   command1="  -dispdev none -e tmp.tcl  > vmd.log 2>&1"
   command1=trim(command1)
   vmd_exe=trim(vmd_exe)
   fcommand=trim(vmd_exe)//trim(command1)
   open(unit=69, file='tmp.sh', status='unknown')
   write(69,*) fcommand
   close(69)
   fcommand="bash tmp.sh"
   call system(fcommand)
   open(unit=109, file='hin_structure.out.cls.color', status='unknown')
   open(unit=110, file='hin_structure.out.cls.lambda', status='unknown')
   write(110,*) "# Time [ps] | N. of mol. into the biggest ice-like cluster"

end subroutine clusters_alloc

subroutine clusters(outxtc,plumed_exe,ns,ws,list_ws,NATOMS,STEP,time, &
                    box_trans,pos,prec,cart,pmpi,cls_stat,r_color,natformat,nat,n_cls_AVE)

   use, intrinsic :: iso_c_binding, only: C_NULL_CHAR, C_PTR, c_f_pointer
   use xtc

   implicit none

   ! Local
   integer :: i, k, STAT_OUT, n_cls, endf, iostat
   integer, allocatable :: a_cls(:)
   type(C_PTR) :: xd_c_out
   type(xdrfile), pointer :: xd_out
   character*100 :: xtcOfile, command1, command2, command3
   character*500 :: fcommand


   ! Arguments
   integer :: cart, pmpi, nat
   integer :: ns, NATOMS, STEP 
   integer, allocatable :: list_ws(:,:), r_color(:)
   real :: time, box_trans(cart,cart), prec, n_cls_AVE
   real, allocatable :: pos(:,:)
   character*4, allocatable :: ws(:)
   character*3 :: outxtc, cls_stat
   character*100 :: plumed_exe, command, natformat


   ! Shitloads of parameters can be tuned when dealing with plumed. However, in here we keep the template fixed,
   ! changing only the indexes of oxygens and hydrogens. If you want to modify plumed behavior, edit plumed.dat_TEMPLATE

   ! Write down the xtc frame of interest for plumed
   xtcOfile='plumed.xtc'
   xtcOfile=trim(xtcOfile)//C_NULL_CHAR
   xd_c_out = xdrfile_open(xtcOfile,"w")
   call c_f_pointer(xd_c_out,xd_out)
   STAT_OUT=write_xtc(xd_out,NATOMS,STEP,time,box_trans,pos,prec)
   STAT_OUT=xdrfile_close(xd_out)

   ! Call plumed executable
   command1="mpirun -np "
   command1=trim(command1)
   write(command2,*) pmpi
   command2=trim(adjustl(command2))
   command3=" --pdb conf.pdb --mf_xtc plumed.xtc > plumed.log 2>&1"
   command3=trim(command3)
   fcommand=trim(command1)//' '//trim(command2)//' '//trim(plumed_exe)//trim(command3)
   call system(fcommand)
   fcommand="tail -1 dfs_surf.dat | cut -d"":"" -f2 | tr ' ' '\n' | sed '/^$/d'  > tmp.dat"
   call system(fcommand)
   open(unit=69, file='tmp.dat', status='old')
   n_cls=0
   endf=0
   do
   read(69,*,iostat=endf)
     if (endf==-1) exit
     n_cls=n_cls+1
   enddo
   allocate(a_cls(n_cls))
   rewind(69)
   do i=1,n_cls
      read(69,*) a_cls(i) ! VMD indexes ! I guess... check ?!
      ! Color ON TOP OF THE RINGS. In hin_structure.out.rings.color you have colors for the rings,
      ! in hin_structure.out.cls.color you have colors for the rings AND for the clusters
      r_color(a_cls(i)+1)=30
   enddo
   close(69)

   ! Analysis of the cluster
   if (trim(adjustl(cls_stat)).eq.'yes') then ! Get asphericity, z-plot...
   endif

   ! Write down the colors for VMD
   write(109,"("//adjustl(natformat)//"i10)") (r_color(k), k=1,nat)

   ! Writeout
   write(110,'(1e10.4,1i10)') time, n_cls

   ! Update average quantities for the cluster
   n_cls_AVE=n_cls_AVE+real(n_cls)

   deallocate(a_cls)

   return

end subroutine clusters

end module MOD_clusters

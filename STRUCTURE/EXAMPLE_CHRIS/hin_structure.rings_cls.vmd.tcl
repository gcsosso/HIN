# Frame number N in this VMD medley corresponds to line number (N+1)*4 in dfs_surf.dat

set sfile md.gro
set tfile hin_structure.out.xtc
# RINGS
#set cfile hin_structure.out.rings.color
# RINGS + PATCH
#set cfile hin_structure.out.rings.color.patch
# RINGS + CLS
set cfile hin_structure.out.cls.color

package require topotools 
package require psfgen 
package require pbctools 

mol load gro $sfile
mol addfile $tfile type xtc waitfor all 0
animate delete beg 0 end 0 
mol delrep 0 top
mol selection {all}
mol material BrushedMetal 
mol addrep top
set topmol [molinfo top]
mol top $topmol
unset topmol
set molid 0
set n [molinfo $molid get numframes]
puts "Reading colors..."
set fp [open $cfile r]
for {set i 0} {$i < $n} {incr i} {
    set chrg($i) [gets $fp]
}
close $fp

proc do_charge {args} {
   global chrg molid
   set f [molinfo $molid get frame]
   set s [atomselect 0 "all"]
   $s set user $chrg($f)
   $s delete
}

trace variable vmd_frame($molid) w do_charge

animate goto start
do_charge

#set nof [molinfo top get numframes]
#puts $nof
#
#set where [pwd]
#
#cd $where
#
#exec rm -r -f snap*
#
#display projection Orthographic
#
#for {set i 0} {$i <= $nof} {incr i} {
#    animate goto $i
#    set idx [ exec echo $i | awk {{printf ("%04d\n",$1)}} ]
#
#    mol delrep all all
#    mol delrep all all
#
#    mol representation CPK 2.000000 0.000000 50.000000 50.000000
#    mol material BrushedMetal
#    mol addrep 0
#    mol modselect 0 0 user > 0
#    mol modcolor 0 0 ColorID 7
#
#    mol representation DynamicBonds 3.200000 0.200000 50.000000
#    mol material BrushedMetal
#    mol addrep 0
#    mol modselect 1 0 user > 0
#    mol modcolor 1 0 ColorID 3
#
#    pbc wrap -center com -centersel "user > 0" 
#    pbc wrap -center com -centersel "user > 0"
#    pbc wrap -center com -centersel "user > 0"
#    pbc wrap -center com -centersel "user > 0"
#
#    display resetview
#
#}

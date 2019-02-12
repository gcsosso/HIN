set sfile md.gro
set tfile hin_structure.out.xtc

# F3
set cfile hin_structure.out.clathrates.f3.color
# F4
#set cfile hin_structure.out.clathrates.f4.color


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
display projection Orthographic
mol modstyle 0 0 CPK
mol modmaterial 0 0 BrushedMetal
mol modselect 0 0 "name OW"
mol modcolor 0 0 User

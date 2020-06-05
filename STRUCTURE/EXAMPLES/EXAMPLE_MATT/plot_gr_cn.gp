set term pos enhanced col dashed eps
set encoding iso_8859_1
set size 0.7,0.7

set outp "atomX.eps"

set xlabel "r (nm)"
#set xtics  20,10,80

set ylabel "g(r)" # tc rgb "#191970"
#set yrange [0:4]
set ytics nomirror

set y2label "Coordination number" rotate by 270
set y2range [0:50]
set y2tics nomirror

#set format y "%2.0t{/ Symbol \327}10^{%L}"
set key off
#set key title "R_m_i_n scaler"
#set key outside 
#set logscale y # ; set logscale y2

plot './hin_structure.out.cryo' u 1:2 title "i001_OO_gr_rCN" w l lt 14 lw 3 lc rgb "#191970" axes x1y1, './hin_structure.out.cryo' u 1:3 w l lt 14 lw 4 lc rgb "#3EB286" axes x1y1, './hin_structure.out.cryo' u 1:4 w l lt 14 lw 4 lc rgb "#FFA500" axes x1y2


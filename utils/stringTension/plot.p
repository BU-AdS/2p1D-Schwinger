#set terminal dumb 120 40
set terminal postscript enhanced
set output '__TITLE__.eps'

set fit quiet
set fit logfile '/dev/null'

set xlabel "t" font ",24"
set ylabel "M_{eff}" font ",24"
set xrange [0:__T__/2+2]
set yrange [0:1.0]

M=1

func1(x)=M

fit [7:15] func1(x) "__MFILE__" using 1:4:5 yerrors via M

set print "chiralExtrap.dat" append  

print __MASS__, M, M_err, FIT_STDFIT*FIT_STDFIT

set label sprintf("mass=%.6f +/- %.6f\n{/Symbol c}^2=%.6f", M, M_err, FIT_STDFIT*FIT_STDFIT) at 3, 0.05

plot [0:__T__/2+1] "__MFILE__" using 1:4:5 with errorbars lc 1 pt 2 ps 3 t "M_{eff}", func1(x) t "effective mass plateau"

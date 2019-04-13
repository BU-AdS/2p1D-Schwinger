#set terminal dumb 120 40
set terminal postscript enhanced
set output 'chiralExtrap.eps'

set xlabel "m_{0}" font ",24"
set ylabel "M_{/Symbol p}" font ",24"
set xrange [-0.2:0.2]
set yrange [0.0:0.8]

A=2.008
B=0.1
C=1.0/3.0
g=0.5
func(x)= A*(((x+B)*(x+B)*(g))**(C))

fit [-0.15:0.05] func(x) "chiralExtrap.dat" using 1:2:3 yerrors via B,g

set label sprintf("M_{/Symbol p} = 2.008*((m_0 + m_{critical})^2g)^{1/3}\n\nm_{critical} = %.6f +/- %.6f\n\ng=%.6f +/- %.6f\n\n{/Symbol b}_{effective} = %.6f\n\n{/Symbol c}^2=%.6f", B, B_err, g, g_err, 1/(g*g), FIT_STDFIT*FIT_STDFIT) at 0.05, 0.3

plot [-0.13:0.20] "chiralExtrap.dat" using 1:2:3 with errorbars lc 1 pt 2 ps 3 t "{/Symbol b}=__BETA__.0 {/Symbol b}_3=1.0", func(x) t "Continuum formula"

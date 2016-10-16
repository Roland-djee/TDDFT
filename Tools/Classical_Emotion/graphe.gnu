reset

set terminal postscript enhanced color
#set terminal postscript color
#set terminal epslatex
#set terminal postscript
set output "Paths.ps"
#set term x11

#set size square

set style line 1 lt 1 lw 1
set style line 2 lt 3 lw 3
set style line 3 lt 4 lw 3

set title 'Classical Electronic trajectories {/Symbol D}t=0,1 [a.u.] {/Symbol w}=0,057 [a.u.] F_0=0,053376 [a.u.]'

#set logscale y

set clip points

set xrange [0:30]
set yrange [0:110]

set xlabel 't_i [a.u.]'
set ylabel 't [a.u.]'

set ticslevel 0

set border 1023-128

set dgrid3d 100,100
#set pm3d
#set hidden3d

set view 60,60

set multiplot

set label 2 at -5,-5,30 center rotate 'Length [a.u.]'
set zrange [0:60]
#set ztics border 1
set grid ztics

splot 'Remoteness.dat' u 1:2:3 notitle w l ls 1

set nolabel
set nogrid

set nolabel

#set ztics nomirror
#set grid ztics
#set label 8 rotate 'Kinetic Energy [a.u.]'
set zrange [0:4]

splot 'Kinetic_Energy.dat' u 1:2:3 notitle w l ls 2

set nolabel
set notics

#set border 1023-128
set zrange [-0.1:0.1]

splot 'Field.dat' u 1:2:3 notitle w l ls 3

set nomultiplot

set out



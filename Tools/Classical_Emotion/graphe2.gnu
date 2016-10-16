reset

set terminal postscript enhanced color
#set terminal postscript color
#set terminal epslatex
#set terminal postscript
set output "Paths.ps"
#set term x11

#set size square

set style line 1 lt 1 lw 1
set style line 2 lt 2 lw 3
set style line 3 lt 3 lw 3

set title 'Classical Electronic trajectories {/Symbol D}t=0,1 [a.u.] {/Symbol w}=0,057 [a.u.] F_0=0,053376 [a.u.]'

set ticslevel 0

set xrange [0:30]
set yrange [0:110]

set xlabel 't_i [a.u.]'
set ylabel 't [a.u.]'

#set border 1023-128

#set pm3d
#set hidden3d

set view 60,60

set multiplot

set label 2 at -5,-5,30 center rotate 'Remoteness [a.u.]'
set zrange [0:60]
set dgrid3d 100,100
set border 1023-128
set grid ztics
set key 50,70,120

splot 'Remoteness.dat' u 1:2:3 title 'Remoteness' w l ls 1

set nolabel
set nogrid
set nodgrid3d
set noborder
set nokey

set label 2 at 32,112,2 center rotate 'Normalized Kinetic Energy  E_k/U_p'
#set origin 0.5,0
set ztics in mirror 0.66666
set zrange [0:4]
set key 50,70,7.7

splot 'Kinetic_Energy.dat' u 1:2:3 title 'Normalized Kinetic Energy' w l ls 2

#set nolabel
#set notics

#set border 1023-128
set zrange [-0.1:0.1]
set key 50,70,0.27

splot 'Field.dat' u 1:2:3 title 'Electric Field' w l ls 3

set nomultiplot

set out



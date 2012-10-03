#!/usr/bin/gnuplot -persist
set terminal png 
set output 'cube-composesquare.png'
set title "Timings for Cubing vs Compose+Square" 
set xlabel "Bits in Discriminant" 
set xrange [14.0000:142.000]
set ylabel "Nanoseconds" 
plot "cube.dat" with lines title 'cube', "compose_square.dat" with lines title 'comp+sqr'
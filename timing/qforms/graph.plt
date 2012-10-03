#!/usr/bin/gnuplot -persist

set terminal png 
set output 's64.png'
set title "Timings for s64 operations" 
set xlabel "Bits in Discriminant" 
set xrange [16.0000:59.000]
set ylabel "Microseconds" 
plot "compose.dat" with lines title 'compose', \
	"square.dat" with lines title 'square', \
	"cube.dat" with lines title 'cube'

set terminal png 
set output 's128.png'
set title "Timings for s128 operations" 
set xlabel "Bits in Discriminant" 
set xrange [60.0000:119.000]
set ylabel "Microseconds" 
plot "compose.dat" with lines title 'compose', \
	"square.dat" with lines title 'square', \
	"cube.dat" with lines title 'cube'
	
set terminal png 
set output 'mpz.png'
set title "Timings for mpz operations" 
set xlabel "Bits in Discriminant" 
set xrange [120.0000:140.000]
set ylabel "Microseconds" 
plot "compose.dat" with lines title 'compose', \
	"square.dat" with lines title 'square', \
	"cube.dat" with lines title 'cube'

set terminal png 
set output 'vs.png'
set title "Timings for compose+square vs cube" 
set xlabel "Bits in Discriminant" 
set xrange [16.0000:140.000]
set ylabel "Microseconds" 
plot "compose_square.dat" with lines title 'comp+sqr', \
	"cube.dat" with lines title 'cube'


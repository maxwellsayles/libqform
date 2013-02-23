#!/usr/bin/gnuplot -persist

set terminal eps
set key left

set output 'compose-all.eps'
set xlabel "Bits in Discriminant" 
set ylabel "Nanoseconds"
plot "compose-64.dat" with lines title '64-bit Multiplication', \
     "compose-128.dat" with lines title "128-bit Multiplication", \
     "compose-mpz.dat" with lines title "GMP Multiplication"

set output 'square-all.eps'
plot "square-64.dat" with lines title '64-bit Squaring', \
     "square-128.dat" with lines title "128-bit Squaring", \
     "square-mpz.dat" with lines title "GMP Squaring"

set output 'cube-all.eps'
plot "cube-64.dat" with lines title '64-bit Cubing', \
     "cube-128.dat" with lines title "128-bit Cubing", \
     "cube-mpz.dat" with lines title "GMP Cubing"


# Cube vs compose+square
set output 'cube-vs-64.eps'
plot "compose_square-64.dat" with lines title '64-bit Multiply w/ Square', \
     "cube-64.dat" with lines title '64-bit Cubing'

set output 'cube-vs-128.eps'
set xrange [59:*]
plot "compose_square-128.dat" with lines title '128-bit Multiply w/ Square', \
     "cube-128.dat" with lines title '128-bit Cubing'
set xrange [*:*]

set output 'cube-vs-mpz.eps'
set xrange [119:*]
plot "compose_square-mpz.dat" with lines title 'GMP Multiply w/ Square', \
     "cube-mpz.dat" with lines title 'GMP Cubing'
set xrange [*:*]
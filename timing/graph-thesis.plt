#!/usr/bin/gnuplot -persist

set terminal eps enhanced
set key left

set output 'compose-all.eps'
set xlabel "Bits in Discriminant" 
set ylabel "Nanoseconds"
set xrange [16:*]
plot "compose-64.dat" with lines title '64-bit Multiplication', \
     "compose-128.dat" with lines title "128-bit Multiplication", \
     "compose-mpz.dat" with lines title "GMP Multiplication", \
     "compose-mpir.dat" with lines title "MPIR Multiplication", \
     "compose-pari.dat" with lines title "PARI Multiplication"

set output 'square-all.eps'
plot "square-64.dat" with lines title '64-bit Squaring', \
     "square-128.dat" with lines title "128-bit Squaring", \
     "square-mpz.dat" with lines title "GMP Squaring", \
     "square-mpir.dat" with lines title "MPIR Squaring", \
     "square-pari.dat" with lines title "PARI Squaring"

set output 'cube-all.eps'
plot "cube-64.dat" with lines title '64-bit Cubing', \
     "cube-128.dat" with lines title "128-bit Cubing", \
     "cube-mpz.dat" with lines title "GMP Cubing", \
     "cube-mpir.dat" with lines title "MPIR Cubing", \
     "cube-pari.dat" with lines title "PARI Cubing"


# Cube vs compose+square
set output 'cube-vs-64.eps'
set xrange [16:59]
plot "compose_square-64.dat" with lines title '64-bit Multiply w/ Square', \
     "cube-64.dat" with lines title '64-bit Cubing'

set output 'cube-vs-128.eps'
set xrange [59:118]
plot "compose_square-128.dat" with lines title '128-bit Multiply w/ Square', \
     "cube-128.dat" with lines title '128-bit Cubing'
set xrange [*:*]

set output 'cube-vs-mpz.eps'
set xrange [119:*]
plot "compose_square-mpz.dat" with lines title 'GMP Multiply w/ Square', \
     "cube-mpz.dat" with lines title 'GMP Cubing'
set xrange [*:*]

set output 'cube-vs-mpir.eps'
set xrange [119:*]
plot "compose_square-mpir.dat" with lines title 'MPIR Multiply w/ Square', \
     "cube-mpir.dat" with lines title 'MPIR Cubing'
set xrange [*:*]

set output 'cube-vs-pari.eps'
set xrange [119:*]
plot "compose_square-pari.dat" with lines title 'PARI Multiply w/ Square', \
     "cube-pari.dat" with lines title 'PARI Cubing'
set xrange [*:*]


# Termination bound.
set xrange [16:59]
set output 'compose-sqrt-vs-64.eps'
plot "compose-sqrt-64.dat" with lines title "Full Square Root", \
     "compose-64.dat" with lines title "Approximate Square Root"
set output 'cube-sqrt-vs-64.eps'
plot "cube-sqrt-64.dat" with lines title "Full Square Root", \
     "cube-64.dat" with lines title "Approximate Square Root"

set xrange [60:118]
set output 'compose-sqrt-vs-128.eps'
plot "compose-sqrt-128.dat" with lines title "Full Square Root", \
     "compose-128.dat" with lines title "Approximate Square Root"
set output 'cube-sqrt-vs-128.eps'
plot "cube-sqrt-128.dat" with lines title "Full Square Root", \
     "cube-128.dat" with lines title "Approximate Square Root"


# Smart 128 bit cubing.
set xrange [60:118]
set output 'smart-cubing-128.eps'
plot 'cube-128.dat' with lines title 'GMP when divisor ≥ 2^{64}', \
     'cube-mixed-mode-arith-128.dat' with lines title 'GMP when product ≥ 2^{128}', \
     'cube-mixed-mode-smart-mod-128.dat' with lines title 'mod when product ≥ divisor'
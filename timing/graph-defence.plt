#!/usr/bin/gnuplot -persist

set terminal eps enhanced
set key left
set xlabel "Bits in Discriminant" 
set ylabel "Average Time (Nanoseconds)"

# 64-bit implementation
set xrange [16:59]
set output 'compose-64.eps'
plot "compose-pari.dat" with lines title "PARI Multiplication", \
     "compose-mpz.dat" with lines title "GMP Multiplication", \
     "compose-64.dat" with lines title '64-bit Multiplication'
set output 'square-64.eps'
plot "square-pari.dat" with lines title "PARI Squaring", \
     "square-mpz.dat" with lines title "GMP Squaring", \
     "square-64.dat" with lines title '64-bit Squaring'
set output 'cube-64.eps'
plot "cube-pari.dat" with lines title "PARI Cubing", \
     "cube-mpz.dat" with lines title "GMP Cubing", \
     "cube-64.dat" with lines title '64-bit Cubing'

# 128-bit implementation
set xrange [60:118]
set output 'compose-128.eps'
plot "compose-pari.dat" with lines title "PARI Multiplication", \
     "compose-mpz.dat" with lines title "GMP Multiplication", \
     "compose-128.dat" with lines title '128-bit Multiplication'
set output 'square-128.eps'
plot "square-pari.dat" with lines title "PARI Squaring", \
     "square-mpz.dat" with lines title "GMP Squaring", \
     "square-128.dat" with lines title '128-bit Squaring'
set output 'cube-128.eps'
plot "cube-pari.dat" with lines title "PARI Cubing", \
     "cube-mpz.dat" with lines title "GMP Cubing", \
     "cube-dyn-128.dat" with lines title '128-bit Cubing'

# Square root approximation.
set xrange [16:59]
set output 'compose-sqrtopt-64.eps'
plot 'compose-sqrt-64.dat' with lines title 'Multiply (Full Sqrt)', \
     'compose-64.dat' with lines title 'Multiply (Approx Sqrt)'
set output 'cube-sqrtopt-64.eps'
plot 'cube-sqrt-64.dat' with lines title 'Cube (Full Sqrt)', \
     'cube-64.dat' with lines title 'Cube (Approx Sqrt)'
     

set xrange [60:118]
set output 'compose-sqrtopt-128.eps'
plot 'compose-sqrt-128.dat' with lines title 'Multiply (Full Sqrt)', \
     'compose-128.dat' with lines title 'Multiply (Approx Sqrt)'
set output 'cube-sqrtopt-128.eps'
plot 'cube-sqrt-128.dat' with lines title 'Cube (Full Sqrt)', \
     'cube-128.dat' with lines title 'Cube (Approx Sqrt)'
     

set xrange [60:118]
set output 'smart-cubing-128.eps'
plot 'cube-128.dat' with lines title 'GMP when divisor ≥ 2^{64}', \
     'cube-mixed-mode-arith-128.dat' with lines title 'GMP when prod ≥ 2^{128}', \
     'cube-mixed-mode-smart-mod-128.dat' with lines title 'mod when prod ≥ divisor'
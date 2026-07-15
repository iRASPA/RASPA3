# Plot the phonon density of states and write a PDF.
# Usage:  gnuplot plot_phonon_dos.gp     -> produces phonon_dos.pdf
# The DOS file has two columns: column 1 = frequency (cm^-1), column 2 = DOS (states per cm^-1 per cell,
# normalized so the integral equals the number of branches 3N). It is integrated over a uniform q-mesh of
# the primitive Brillouin zone, so it is directly comparable to an inelastic-neutron generalized phonon DOS.

file = 'output/phonon_dos.s0.txt'

set terminal pdfcairo enhanced font 'Helvetica,18' size 8in,6in linewidth 2
set output 'phonon_dos.pdf'

set xlabel 'Frequency (cm^{-1})' font 'Helvetica,20'
set ylabel 'Phonon DOS (states / cm^{-1})' font 'Helvetica,20'
set title  'IRMOF-1 phonon density of states (primitive cell)' font 'Helvetica,22'
set tics   font 'Helvetica,16'
set grid ytics lt 1 lc rgb 'gray80'
set xzeroaxis lt 1 lc rgb 'gray50'
set xrange [*:*]
set yrange [0:*]
unset key

set style fill solid 0.35 noborder
plot file using 1:2 with filledcurves y1=0 lc rgb '#4060c0', \
     file using 1:2 with lines lw 2 lc rgb 'navy'

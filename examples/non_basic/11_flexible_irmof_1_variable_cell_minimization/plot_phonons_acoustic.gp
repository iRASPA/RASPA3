# Plot the acoustic phonon branches (the lowest modes that vanish at Gamma) and write a PDF.
# Usage:  gnuplot plot_phonons_acoustic.gp   -> produces phonon_dispersion_acoustic.pdf
# The band file stores frequencies sorted ascending per k-point: column 1 is the path
# coordinate, columns 2..N+1 are the bands, and high-symmetry rows end in "# LABEL".

file = 'output/phonon_dispersion.s0.txt'

# --- larger fonts throughout; pdfcairo gives crisp vector output ---
set terminal pdfcairo enhanced font 'Helvetica,18' size 8in,6in linewidth 2
set output 'phonon_dispersion_acoustic.pdf'

# high-symmetry labels -> xtics (data rows that carry a trailing "# LABEL")
tics = system("awk '!/^#/ && /#/ {printf \"%s\\\"%s\\\" %s\", sep, $NF, $1; sep=\", \"}' " . file)
eval 'set xtics (' . tics . ')'

set grid xtics lt 1 lc rgb 'gray70'
set xzeroaxis lt 1 lc rgb 'gray50'      # zero line: acoustic branches touch it at Gamma
set ylabel 'Frequency (cm^{-1})' font 'Helvetica,20'
set title  'IRMOF-1 acoustic phonon branches' font 'Helvetica,22'
set tics   font 'Helvetica,16'
set xrange [0:*]
unset key

# columns 2:4 = the three lowest (acoustic) branches that go to zero at Gamma
plot for [i=2:4] file using 1:i with lines lw 3

# Plot the full phonon dispersion (all branches) and write a PDF.
# Usage:  gnuplot plot_phonons.gp        -> produces phonon_dispersion.pdf
# The band file stores frequencies sorted ascending per k-point: column 1 is the path
# coordinate, columns 2..N+1 are the bands, and high-symmetry rows end in "# LABEL".
# The dispersion is computed on the primitive cell (true unfolded band structure), so with
# 106 atoms in the IRMOF-1 primitive cell there are 3*106 = 318 branches.

file = 'output/phonon_dispersion.s0.txt'

# --- larger fonts throughout; pdfcairo gives crisp vector output ---
set terminal pdfcairo enhanced font 'Helvetica,18' size 8in,6in linewidth 2
set output 'phonon_dispersion.pdf'

# number of band columns (use an interior row with no trailing '# LABEL' annotation)
nbands = int(system("awk '!/#/ && NF>1 {print NF-1; exit}' " . file))

# high-symmetry labels -> xtics (data rows that carry a trailing "# LABEL")
tics = system("awk '!/^#/ && /#/ {printf \"%s\\\"%s\\\" %s\", sep, $NF, $1; sep=\", \"}' " . file)
eval 'set xtics (' . tics . ')'

set grid xtics lt 1 lc rgb 'gray70'
set xzeroaxis lt 1 lc rgb 'gray50'      # zero line: acoustic branches touch it at Gamma
set ylabel 'Frequency (cm^{-1})' font 'Helvetica,20'
set title  'IRMOF-1 phonon dispersion (primitive cell)' font 'Helvetica,22'
set tics   font 'Helvetica,16'
set xrange [0:*]
unset key

# all branches (columns 2 .. nbands+1)
plot for [i=2:nbands+1] file using 1:i with lines lw 1 lc rgb 'navy'

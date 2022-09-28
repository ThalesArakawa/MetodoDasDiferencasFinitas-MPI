set terminal png enhanced font "utf8"
set title "SpeedUp MPI nó único"
set xlabel "N. threads"
set ylabel "Eficiência"
set xrange [0:49]
set key right center
set output "BenchmarkSiestaMPIEficiência.png"
plot "speedUp107.dat" using 1:2 title "I7-2600" with linespoints pt 6 linecolor 150, "speedUp1no710.dat" using 1:2 title "E5-2695v2" with linespoints pt 6 linecolor 100, "speedUp1noSeq.dat" using 1:2 title "CLG 6252" with linespoints pt 6 linecolor 10

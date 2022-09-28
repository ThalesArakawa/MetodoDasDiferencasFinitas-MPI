set terminal png enhanced font "utf8"
set title "Comparativo de Processadores MPI/Nó único"
set xlabel "Número de Trabalhadores"
set ylabel "Tempo de Execução (s)"
set xrange [0:49]
set key right bottom
set output "comparatiboOmpMpi.png"
plot "tempo107MPI.txt" using 1:2 title "I7-2600" with linespoints pt 1, "tempo1no710.dat" using 1:2 title "E5-2695v2" with linespoints pt 2 , "tempo1noSeq.dat" using 1:2 title "CLG 6252" with linespoints pt 3

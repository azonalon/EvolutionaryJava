#!/bin/bash
javac -g physics/*.java &&

java -ea physics.SoftBody "beam oscillation rotation" 0.01 100  > tempfile.dat &&
gnuplot -e "set parametric;\
set terminal qt 0;\
plot 'tempfile.dat' u 7:(column(8)) w linespoints ls 0;\
set terminal qt 1;\
plot for [n=7:8] 'tempfile.dat' u 1:(column(n)) w linespoints ls n;\
set terminal qt 2;\
plot for [n=10:11] 'tempfile.dat' u 1:(column(n)) w linespoints ls n;\
set terminal qt 3;\
plot 'tempfile.dat' u 1:(column(12)) w linespoints ls 1;\
pause -1"

if false; then
    java -ea physics.SoftBody "orbit" 0.1 800 > tempfile.dat &&
    gnuplot -e "set parametric;\
    plot for [n=2:4:2] 'tempfile.dat' u n:(column(n+1)) ls n pt 1  ps 1;\
    set terminal qt 1;\
    plot for [n=6:6:1] 'tempfile.dat' u 1:(column(n)) ls n pt 1  ps 1;\
    pause -1"

    java -ea physics.SoftBody "beam oscillation" 0.1 500  > tempfile.dat &&
    gnuplot -e "set parametric;\
    set terminal qt 0;\
    plot for [n=5:7] 'tempfile.dat' u 1:(column(n)) w linespoints ls n  title 'y th l';\
    set terminal qt 1;\
    plot for [n=8:8] 'tempfile.dat' u 1:(column(n)) w linespoints ls n  title 'energy';\
    pause -1"



    java -ea physics.SoftBody "basic rotation" 0.01 400 > tempfile.dat &&
    gnuplot -e "set parametric;\
    plot for [n=2:4:2] 'tempfile.dat' u n:(column(n+1)) ls n pt 1  ps 1;\
    set terminal qt 1;\
    plot for [n=6:6:1] 'tempfile.dat' u 1:(column(n)) ls n pt 1  ps 1;\
    pause -1"

    java -ea physics.SoftBody "relative motion damping" 0.01 1000 > tempfile.dat &&
    gnuplot -e "\
    plot for [n=2:3] 'tempfile.dat' u 1:(column(n)) ls n pt 1  ps 1;\
    set terminal qt 1;\
    plot for [n=4:4:1] 'tempfile.dat' u 1:(column(n)) ls n pt 1  ps 1;\
    pause -1"

    java -ea physics.SoftBody "relative motion" 0.1 100 > tempfile.dat &&
    gnuplot -e "\
    plot for [n=2:3] 'tempfile.dat' u 1:(column(n)) ls n pt 1  ps 1;\
    set terminal qt 1;\
    plot for [n=4:4:1] 'tempfile.dat' u 1:(column(n)) ls n pt 1  ps 1;\
    pause -1"
fi
rm -f tempfile

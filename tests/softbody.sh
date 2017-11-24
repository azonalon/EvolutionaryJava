#!/bin/bash
javac -g physics/*.java &&
javac -g physics/*/*.java

if [ 0 -eq $? ]; then
    java -ea physics.test.TwoCellTests "damped symmetric rotation oscillation" 0.1 1000  > tempfile.dat
    gnuplot -e "set parametric;\
    set terminal qt 0;\
    plot 'tempfile.dat' u 7:(column(8)) w linespoints ls 1,\
         'tempfile.dat' u 2:(column(3)) w linespoints ls 2;\
    set terminal qt 1;\
    plot 'tempfile.dat' u 1:(column(4) - column(13)) w linespoints ls 3;\
    pause -1"

    exit 0

    java -ea physics.test.TwoCellTests "symmetric rotation oscillation" 0.01 840  > tempfile.dat
    gnuplot -e "set parametric;\
    set terminal qt 0;\
    plot 'tempfile.dat' u 7:(column(8)) w linespoints ls 1, \
         'tempfile.dat' u 2:(column(3)) w linespoints ls 2 ;\
    set terminal qt 1;\
    plot 'tempfile.dat' u 1:(column(12)) w linespoints ls 0, \
         'tempfile.dat' u 1:(2.4*(column(13) - column(4))**2 - 0.035) w linespoints ls 1; \
    set terminal qt 2;\
    plot 'tempfile.dat' u 1:((column(13) - column(4))) w linespoints ls 0; \
    pause -1"


    java -ea physics.test.TwoCellTests "beam oscillation" 0.1 500  > tempfile.dat &&
    gnuplot -e "set parametric;\
    set terminal qt 0;\
    plot for [n=5:7] 'tempfile.dat' u 1:(column(n)) w linespoints ls n  title 'y th l';\
    set terminal qt 1;\
    plot for [n=8:8] 'tempfile.dat' u 1:(column(n)) w linespoints ls n  title 'energy';\
    pause -1"

    java -ea physics.test.TwoCellTests "beam oscillation rotation" 0.001 1000  > tempfile.dat &&
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

    java -ea physics.test.TwoCellTests "orbit" 0.1 800 > tempfile.dat &&
    gnuplot -e "set parametric;\
    plot for [n=2:4:2] 'tempfile.dat' u n:(column(n+1)) ls n pt 1  ps 1;\
    set terminal qt 1;\
    plot for [n=6:6:1] 'tempfile.dat' u 1:(column(n)) ls n pt 1  ps 1;\
    pause -1"




    java -ea physics.test.TwoCellTests "basic rotation" 0.01 400 > tempfile.dat &&
    gnuplot -e "set parametric;\
    plot for [n=2:4:2] 'tempfile.dat' u n:(column(n+1)) ls n pt 1  ps 1;\
    set terminal qt 1;\
    plot for [n=6:6:1] 'tempfile.dat' u 1:(column(n)) ls n pt 1  ps 1;\
    pause -1"

    java -ea physics.test.TwoCellTests "relative motion damping" 0.01 1000 > tempfile.dat &&
    gnuplot -e "\
    plot for [n=2:3] 'tempfile.dat' u 1:(column(n)) ls n pt 1  ps 1;\
    set terminal qt 1;\
    plot for [n=4:4:1] 'tempfile.dat' u 1:(column(n)) ls n pt 1  ps 1;\
    pause -1"

    java -ea physics.test.TwoCellTests "relative motion" 0.1 100 > tempfile.dat &&
    gnuplot -e "\
    plot for [n=2:3] 'tempfile.dat' u 1:(column(n)) ls n pt 1  ps 1;\
    set terminal qt 1;\
    plot for [n=4:4:1] 'tempfile.dat' u 1:(column(n)) ls n pt 1  ps 1;\
    pause -1"
fi
rm -f tempfile

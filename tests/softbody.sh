#!/bin/bash
javac physics/*.java &&

java -ea physics.SoftBody "beam oscillation" 0.01 200 > tempfile.dat &&
gnuplot -e "set parametric;\
plot for [n=2:4:2] 'tempfile.dat' u n:(column(n+1)) ls n pt 1  ps 1;\
pause -1"

exit 0

java -ea physics.SoftBody "orbit" 0.01 800 > tempfile.dat &&
gnuplot -e "set parametric;\
plot for [n=2:4:2] 'tempfile.dat' u n:(column(n+1)) ls n pt 1  ps 1;\
pause -1"

java -ea physics.SoftBody "basic rotation" 0.01 400 > tempfile.dat &&
gnuplot -e "set parametric;\
plot for [n=2:4:2] 'tempfile.dat' u n:(column(n+1)) ls n pt 1  ps 1;\
pause -1"

java -ea physics.SoftBody "relative motion damping" 0.01 1000 > tempfile.dat &&
gnuplot -e "\
plot for [n=2:3] 'tempfile.dat' u 1:(column(n)) ls n pt 1  ps 1;\
pause -1"

java -ea physics.SoftBody "relative motion" 0.1 100 > tempfile.dat &&
gnuplot -e "\
plot for [n=2:3] 'tempfile.dat' u 1:(column(n)) ls n pt 1  ps 1;\
pause -1"

rm -f tempfile

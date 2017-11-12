#!/bin/sh
javac physics/*.java &&
java physics.Collision &&
gnuplot -e "\
plot for [n=2:6] '< java -ea physics.Collision' u 1:(column(n)) ;\
pause -1\
"

#!/bin/bash
javac physics/test/*.java &&
gnuplot -e "\
plot for [n=2:2] '< java -ea physics/test/Leapfrog' u 1:(column(n));\
pause -1\
"

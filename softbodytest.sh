#!/bin/bash
javac physics/*.java &&
java physics.SoftBody
# gnuplot -e "\
# plot for [n=2:6] '< java -ea physics.Collision' u 1:(column(n)) ;\

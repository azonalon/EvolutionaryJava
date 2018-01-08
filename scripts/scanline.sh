#!/bin/sh
rm -rf build/test-results/physics/ScanLine;
gradle cleanTest &&
gradle test --console rich --tests physics.ElasticModelTests.backwardEulerTest &&
gnuplot -e "\
plot 'build/test-results/physics/ScanLine/scanline_${1}_${2}.dat'\
 u 1:3, '' u 1:2, '' u 1:4;\
 pause -1"

#!/bin/sh
# rm -rf build/test-results/physics/ScanLine
# gradle test --console rich --tests physics.ElasticModelTests.backwardEulerTestMultiTriangle
gnuplot -e "\
plot 'build/test-results/physics/ScanLine/scanline_${1}_${2}.dat'\
 u 1:3, '' u 1:2, '' u 1:4, '' u 1:(column(4)-column(3));\
 pause -1"

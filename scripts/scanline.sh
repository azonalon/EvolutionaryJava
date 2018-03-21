#!/bin/sh
# rm -rf build/test-results/physics/ScanLine
# gradle test --console rich --tests physics.ElasticModelTests.backwardEulerTestMultiTriangle
gnuplot -e "\
plot 'build/test-results/physics/ScanLine/scanline_${1}_${2}.dat'\
 u 1:3 title \"alpha\", '' u 1:2 title \"phi\", '' u 1:3 title \"dphi\", '' u 1:4 title \"dphiapprox\";\
 pause -1"

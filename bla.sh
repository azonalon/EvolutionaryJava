#!/bin/sh
# rm -rf build/test-results/physics/SoftBody/ScanLine &&
#     gradle cleanTest &&
#     gradle --console rich test --tests physics.$1;
#     cd build/test-results/physics/SoftBody/ScanLine/ &&
#     gnuplot &&
#     cd -
#
#

rm -rf build/test-results/physics/ElasticModelTests
gradle cleanTest &&
gradle test --console rich --tests physics.ElasticModelTests.backwardEulerTest &&
# cat 'build/test-results/physics/ElasticModelTests/backwardEuler1.dat' ;
gnuplot -e "\
set parametric;\
plot for[i=2:6:2] \
'build/test-results/physics/ElasticModelTests/backwardEuler1.dat'\
 u (column(i)):(column(i+1));\
 pause -1"

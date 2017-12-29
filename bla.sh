#!/bin/sh
rm -rf build/test-results/physics/SoftBody/ScanLine &&
    gradle cleanTest &&
    gradle --console rich test --tests physics.$1;
    cd build/test-results/physics/SoftBody/ScanLine/ &&
    gnuplot &&
    cd -

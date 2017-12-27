#!/bin/sh
javac ArcTan.java && java ArcTan > atan.dat
gnuplot -e 'plot "atan.dat" u 1:2 title "atan(y/x)",
 "atan.dat" u 1:3 title "atan2(y, x)",
 "atan.dat" u 1:4 title "2 atan(y/(1+x))";
pause -1;
'

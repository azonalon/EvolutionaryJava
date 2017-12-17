#!/bin/sh
javac ArcTan.java && java ArcTan > atan.dat
gnuplot -e 'plot "atan.dat" u 1:2 title "atan(x/y)",
 "atan.dat" u 1:3 title "atan2(x, y)",
 "atan.dat" u 1:4 title "atan(y/(1+x))";
pause -1;
'

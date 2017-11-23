#!/bin/python
import sys
import time
import subprocess as sp
import subprocess
import pty
import os

i, o = pty.openpty()

proc = subprocess.Popen(['gnuplot'], shell=True, stdin=i)
os.write(i, b"plot sin(x);pause 4\n")
# print(str(os.read(o, 100)))
while True:
    print('write please')
    t = sys.stdin.readline().encode()
    os.write(i, t)
    print(t, 'written')

import math
print(math.atan2(0, 1));
print(math.atan2(1, 0));
print(math.atan2(-0.01, -1));
print(math.atan2(-1, 0));

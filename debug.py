#!/bin/python3
import sys
import time
import subprocess as sp
import subprocess
import pty
import os
from scipy.optimize import minimize, basinhopping

def f(x):
    r=sp.run(["java", "-ea", "physics.SoftBody", "minimize",
              str(x[0]),str(x[1]),str(x[2]),str(x[3])],
             stdout=sp.PIPE, stderr=sp.PIPE)
    return float(r.stdout.split(b"\n")[-2])

min = minimize(f, x0=[1, 1, 1, 1],
         bounds=((0, 10),(0, 10), (0, 10), (0, 10)) )
min2 = basinhopping(f, x0=[1, 1, 1, 1])
f(min.x)
f([0, 0, 0.2, 0.2])
min.x
f(0, 2, 0, 3)

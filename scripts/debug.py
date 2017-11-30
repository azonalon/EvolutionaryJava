#!/bin/python3
import subprocess as sp
import re
from scipy.optimize import minimize, basinhopping

def f(x):
    r=sp.run(["java", "-ea", "physics.SoftBody", "minimize",
              str(x[0]),str(x[1]), "0", "0"],
             stdout=sp.PIPE, stderr=sp.PIPE)
    return float(r.stdout.split(b"\n")[-2])

min = minimize(f, x0=[2.4, -0.034],
         bounds=[(0, 10),(-1, 1)])
print(min)
min2 = basinhopping(f, x0=[1, 0],
                    callback=lambda x,f,a: print("Minima found", x, f))
f([2.4, 2.4 * -0.035])
f([0, 0])
f(min.x)
min
f([0, 0, 0.2])
min.x
f(0, 2, 0, 3)
[a for a in range(20) if a < 10]

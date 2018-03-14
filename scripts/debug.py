#!/bin/python3
import subprocess as sp
import re
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.optimize import minimize, basinhopping
# S = np.linspace(-np.pi, np.pi, 100)
# X = np.cos(S)
# Y = np.sin(S)
# plt.plot(S, np.arctan2(X, Y))
# plt.plot(S, np.arctan(X/Y))
# plt.show()

from numpy import asarray, diag, sqrt, hypot, array, sin
from math import atan2
from numpy.linalg import svd


def sqrd(x):
    return x*x

def svd2(m):
    y1, x1 = (m[1, 0] + m[0, 1]), (m[0, 0] - m[1, 1])
    y2, x2 = (m[1, 0] - m[0, 1]), (m[0, 0] + m[1, 1])


    h1 = hypot(y1, x1)
    t1 = x1 / h1

    h2 = hypot(y2, x2)
    t2 = x2 / h2

    cc = sqrt((1 + t1) * (1 + t2))
    ss = sqrt((1 - t1) * (1 - t2))
    cs = sqrt((1 + t1) * (1 - t2))
    sc = sqrt((1 - t1) * (1 + t2))

    c1, s1 = (cc - ss) / 2, (sc + cs) / 2,
    u1 = asarray([[c1, -s1], [s1, c1]])

    d = asarray([(h1 + h2) / 2, (h1 - h2) / 2])
    sigma = diag(d)

    if h1 != h2:
        u2 = diag(1 / d).dot(u1.T).dot(m)
    else:
        u2 = diag([1 / d[0], 0]).dot(u1.T).dot(m)

    return u1, d, u2

for _ in range(1):
    x = np.random.random((2,2))
    x = np.diag([3.1,3.1])
    sv2=svd(x)
    sv1=svd2(x)
    print("M:", x, "Det", np.linalg.det(x))
    print("U:", sv1[0], sv2[0], sep='\n')
    print("S:", sv1[1], sv2[1], sep='\n')
    print("V:", sv1[2], sv2[2], sep='\n')
    print(hypot(0,0))

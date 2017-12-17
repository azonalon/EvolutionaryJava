#!/bin/python3
import subprocess as sp
import re
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.optimize import minimize, basinhopping
S = np.linspace(-np.pi, np.pi, 100)
X = np.cos(S)
Y = np.sin(S)
plt.plot(S, np.arctan2(X, Y))
plt.plot(S, np.arctan(X/Y))
plt.show()

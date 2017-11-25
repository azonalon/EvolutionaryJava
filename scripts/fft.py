#!/bin/python
import numpy.fft as fft
from scipy.fftpack import rfft, rfftfreq
import numpy as np
import tempfile, sys

with tempfile.NamedTemporaryFile(mode='w') as f:
    f.write(sys.stdin.read())
    f.flush()
    arr = np.genfromtxt(f.name, delimiter=' ', unpack=True)

n=len(arr[0])
arr[0] = rfftfreq(n, (arr[0,-1] - arr[0,0])/n)
for i in range(1, len(arr)):
    arr[i] = rfft(arr[i])

for line in np.transpose(arr):
    print(" ".join("% 10.4f" % abs(v) for v in line))

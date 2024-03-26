#!/bin/python3


import numpy as np
import scipy.linalg

A = np.random.random((8,8)) + 1j*np.random.random((8,8))
B = np.random.random((8,8)) + 1j*np.random.random((8,8))



for i in range(8):
    for j in range(8):
        print(np.real(A[i,j]),np.imag(A[i,j]))

print("================================================================")

for i in range(8):
    for j in range(8):
        print(np.real(B[i,j]),np.imag(B[i,j]))

print("================================================================")

w,v = scipy.linalg.eig(A,b=B)

for i in range(8):
    print(np.real(w[i]),np.imag(w[i]))
print("================================================================")
for i in range(8):
    for j in range(8):
        print(np.real(v[i,j]),np.imag(v[i,j]))

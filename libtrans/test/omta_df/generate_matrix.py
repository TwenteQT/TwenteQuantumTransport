#!/bin/python3


import numpy as np

A = np.random.random((8,8)) + 1j*np.random.random((8,8))

B = np.linalg.inv(A) 


for i in range(8):
    for j in range(8):
        print(np.real(A[i,j]),np.imag(A[i,j]))

print("================================================================")

for i in range(8):
    for j in range(8):
        print(np.real(B[i,j]),np.imag(B[i,j]))


print(np.dot(A,B))

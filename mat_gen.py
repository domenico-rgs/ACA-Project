#!/usr/bin/env python3
import sys
import numpy as np

def generate_matrix(nrows, ncols):
    return 2.5*np.random.randn(nrows, ncols)

for i in range(1,len(sys.argv)):
    nrows, ncols = sys.argv[i].split("x")
    matrix = generate_matrix(int(nrows), int(ncols))
    with open("mat"+str(i)+".txt",'wb') as f:
        f.write((" ".join([nrows, ncols])+"\n").encode())
        np.savetxt(f, matrix, fmt='%.5f')

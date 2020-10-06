import numpy as np


def tridiagonal(A, f):
    n = len(A)
    a = np.zeros(n)
    b = np.zeros(n - 1)
    c = np.zeros(n - 1)
    d = f.copy()
    a = A.diagonal().copy()
    b = A.diagonal(offset=1).copy()
    c = A.diagonal(offset=-1).copy()
    for i in range(n-1):
        b[i] = b[i] / a[i]
        a[i + 1] = a[i + 1] - c[i] * b[i]
    d[0] = d[0] / a[0]
    for j in range(1, n):
        d[j] = (d[j] - c[j-1] * d[j - 1]) / a[j]
    for k in range(n - 2, -1, -1):
        d[k] = d[k] - b[k] * d[k + 1]

    return d

#def LUbandmatrix(A, r, s):
#    n = A.shape[1]
#    m = r+s+1


A = [[8., -1., 0., 0., 0.],
     [2., 8., 1., 0., 0.],
     [0., 1., 8., -1., 0.],
     [0., 0., 2., 8., 1.],
     [0., 0., 0., 1., 8.]]
A = np.array(A)
b = [2.2, -2.54, 8.26, 8.32, -4.3]
b = np.array(b)
print(tridiagonal(A, b))
LUbandmatrix(A,1,1)

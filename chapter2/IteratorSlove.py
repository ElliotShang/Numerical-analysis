import numpy as np
from numpy.core._multiarray_umath import ndarray


def Jacobi_nonmatrix(A, b):
    n = len(A)
    x = np.zeros(n)
    y = np.ones(n)
    N = 100
    tmp = np.zeros(n)
    #    for k in range(N):
    while max(abs((y - tmp))) > 0.001:
        for i in range(0, n):
            s1 = sum(A[i][t] * x[t] for t in range(i))
            s2 = sum(A[i][t] * x[t] for t in range(i + 1, n))
            y[i] = (b[i] - s1 - s2) / A[i][i]
        tmp = x.copy()
        x = y.copy()  # You can understand x and y as pointers to a specific array.
        # Direct assignment will cause x and y to point to the same object, causing the x and y values to change
        # together during the iteration.
    return y


def Jacobi_matrix(A, b):
    # Use numpy matrix operation, save time
    x = np.zeros(len(b))
    D = np.diagflat(np.diag(A))
    LU = A - D
    N = 100
    for k in range(N):
        x = np.dot(np.dot(-np.linalg.inv(D), LU), x) + np.dot(np.linalg.inv(D), b)
    #    D = np.diag(A)
    #    LU = A - np.diagflat(D)
    #    N = 100
    #    for k in range(N):
    #       x = (b - np.dot(LU, x))/D
    return x


def Gauss_seidel(A, b):
    n = len(A)
    x = np.zeros(n)
    tmp = np.ones(n)
    #   N = 100
    #   for k in range(N):
    while max(abs((x - tmp))) > 0.001:
        tmp = x.copy()
        for i in range(n):
            x[i] = b[i] / A[i][i] - sum((A[i][j] / A[i][i]) * x[j] for j in range(i)) - sum(
                (A[i][j] / A[i][i]) * x[j] for j in range(i + 1, n))
    return x




# Test case 1, problem 15, problem 17

A = [[5., 2., 1.], [-1., 4., 2.], [2., -5., 10.]]
A = np.array(A)
b = [-12., 10., 1.]
b = np.array(b)
print(Jacobi_nonmatrix(A, b))
print(Jacobi_matrix(A, b))
print(Gauss_seidel(A, b))

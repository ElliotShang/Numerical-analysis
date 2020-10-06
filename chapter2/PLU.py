import numpy as np


def LU_partial_decomposition(matrix):
    n, m = matrix.shape
    P = np.identity(n)
    L = np.identity(n)
    U = matrix.copy()
    PF = np.identity(n)
    LF = np.zeros((n, n))
    for k in range(0, n - 1):
        index = np.argmax(abs(U[k:, k]))
        index = index + k
        if index != k:
            P = np.identity(n)
            P[[index, k], k:n] = P[[k, index], k:n]
            U[[index, k], k:n] = U[[k, index], k:n]
            PF = np.dot(P, PF)
            LF = np.dot(P, LF)
        L = np.identity(n)
        for j in range(k + 1, n):
            L[j, k] = -(U[j, k] / U[k, k])
            LF[j, k] = (U[j, k] / U[k, k])
        U = np.dot(L, U)
    np.fill_diagonal(LF, 1)
    return PF, LF, U


A = [[1, 8, 2, 3],
     [-6, -3, 8, 1],
     [2, 4, 4, 2],
     [10, 5, -5, 6]]
A = np.array(A)
P1, L1, U1 = LU_partial_decomposition(A)
print(P1)
print(L1)
print(U1)

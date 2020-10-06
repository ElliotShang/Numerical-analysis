import numpy as np
from scipy.linalg import lu, lu_factor

def Doolittle_withoutprivot(A):
    n = len(A)
    M = np.zeros((n, n))
    #   U = np.zeros((n, n))
    for k in range(n):
        # upper triangular
        s = 0
        for j in range(k, n):
            s = sum(M[k][t] * M[t][j] for t in range(k))
            M[k][j] = (A[k][j] - s)

        # lower triangular
        for i in range(k + 1, n):
            #        if i == k:
            #            A[i][i] = 1
            #       else:
            s = sum(M[i][t] * M[t][k] for t in range(k))
            M[i][k] = ((A[i][k] - s) / M[k][k])

    print(M)
    return M


def Doolittleprivot(A):
    n = len(A)
    lu = np.zeros((n, n))
    M = np.zeros(n)
    for k in range(n):
        s = np.zeros(n)
        max = -1e100
        for i in range(k, n):
            s[i] = A[i][k] - sum(lu[i][t] * lu[t][k] for t in range(k))
        for i in range(n):
            if max < abs(s[i]):
                maxindex = i
                max = abs(s[i])
        for t in range(k):
            tmp = lu[k][t]
            lu[k][t] = lu[maxindex][t]
            lu[maxindex][t] = tmp
#        lu[k], lu[maxindex] = lu[maxindex], lu[k]
#        lu[[k, maxindex]] = lu[[maxindex, k]]
        for t in range(k, n):
            tmp = A[k][t]
            A[k][t] = A[maxindex][t]
            A[maxindex][t] = tmp
#        A[k], A[maxindex] = A[maxindex], A[k]
#        A[[k, maxindex]] = A[[maxindex, k]]
        tmp = s[k]
        s[k] = s[maxindex]
        s[maxindex] = tmp
        lu[k][k] = s[k]
        for j in range(k+1, n):
            lu[k][j] = A[k][j] - sum(lu[k][t]*lu[t][j] for t in range(k))
        for i in range(k+1, n):
            lu[i][k] = s[i]/lu[k][k]

    print(lu)
    return lu


def SolveLU(A, b):
    n = len(b)
    x = np.zeros(n)
    y = np.zeros(n)
    y[0] = b[0]
    for i in range(1, n):
        y[i] = b[i] - sum(A[i][t] * y[t] for t in range(i))
    x[n - 1] = y[n - 1] / A[n - 1][n - 1]
    for i in range(n - 2, -1, -1):
        x[i] = (y[i] - sum(A[i][t] * x[t] for t in range(i + 1, n))) / A[i][i]
    print(x)
    return x


mat = [[1, 8, 2, 3],
       [-6, -3, 8, 1],
       [2, 4, 4, 2],
       [10, 5, -5, 6]]
#b = [2.5, 1.8, 7.2]
print('主元')
Doolittleprivot(mat)
print('无主元')
Doolittle_withoutprivot(mat)
print('标准库')
P,L,U  = lu(mat)
print(P)
print(L)
print(U)
#SolveLU(A, b)

import numpy as np


def GEPP(A, b):
    n = len(A)
    for i in range(n - 1):
        max = -1e100
        for r in range(i, n):
            if max < abs(A[r][i]):
                maxindex = r
                max = abs(A[r][i])
        if maxindex != i:
            A[[i, maxindex]] = A[[maxindex, i]]
            b[[i, maxindex]] = b[[maxindex, i]]
        for j in range(i + 1, n):
            m = A[j][i] / A[i][i]
            A[j][i] = m
            for k in range(i + 1, n):
                A[j][k] = A[j][k] - m * A[i][k]
            b[j] = b[j] - m * b[i]
    print(A)
    print(b)
    x = np.zeros(n)
    k = n - 1
    x[k] = b[k] / A[k][k]
    while k >= 0:
        x[k] = (b[k] - np.dot(A[k][k + 1:], x[k + 1][:])) / A[k][k]
        k = k - 1
    return x


A = np.array([[0.5, 1.1, 3, 1], [5, 0.96, 6.5], [2, 4.5, 0.36]])
b = np.array([[6.], [0.96], [2.]])
print(GEPP(A, b))

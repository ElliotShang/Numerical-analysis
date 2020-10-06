import numpy as np


def gaussian_elimination_with_pivot(A, b):
    n = len(A)
    for k in range(n):
        max = -1e100
        for p in range(k, n):
            if max < abs(A[p][k]):
                max_row = p
                max = abs(A[p][k])
    # A[p], A[max_row] = A[max_row], A[p]
        tmp = A[max_row]
        A[max_row] = A[k]
        A[k] = tmp
        b[k], b[max_row] = b[max_row],b[k]
        print(A)
        print(b)
        for i in range(k + 1, n):
            m = A[i][k] / A[k][k]
            for j in range(k, n):
                A[i][j] = A[i][j] - m * A[k][j]
            b[i] = b[i] - m * b[k]

    print(A)
    # back substitution
    x = [0] * n
   #x[n-1] = b[n-1] / A[n-1][n-1]
    for i in range(n-1, -1, -1):
        s = sum(A[i][j] * x[j] for j in range(i, n))
        x[i] = (b[i] - s) / A[i][i]
    return x


if __name__ == '__main__':
    # m = [[0,-2,6,-10], [-1,3,-6,5], [4,-12,8,12]]
    # m = [[1,-1,3,2], [3,-3,1,-1], [1,1,0,3]]
    m = [[0.5, 1.1, 3.1], [5, 0.96, 6.5], [2, 4.5, 0.36]]
    b = [6, 0.96, 0.02]
    print(gaussian_elimination_with_pivot(m, b))

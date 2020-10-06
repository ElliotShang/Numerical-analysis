import numpy as np


def gaussian_elimination_with_pivot(m):

    # forward elimination
    n = len(m)
    print(n)
    for k in range(n):
        max = -1e100
        for r in range(k, n):
            if max < abs(m[r][k]):
                max_row = r
                max = abs(m[r][k])
        m[k], m[max_row] = m[max_row], m[k]
        print(m)
        for i in range(k+ 1, n):
            m[i] = [m[i][j] - m[k][j] * m[i][k] / m[k][k] for j in range(n + 1)]

    if m[n - 1][n - 1] == 0: raise ValueError('No unique solution')

    print(m)

    # backward substitution
    x = [0] * n
    for i in range(n - 1, -1, -1):
        s = sum(m[i][j] * x[j] for j in range(i, n))
        x[i] = (m[i][n] - s) / m[i][i]
    return x


'''
# shorter way to pivot but cannot run in trinket
def pivot(m, n, i):
  max_row = max(range(i, n), key=lambda r: abs(m[r][i]))
  m[i], m[max_row] = m[max_row], m[i]
'''


def pivot(m, n, i):
    max = -1e100
    for r in range(i, n):
        if max < abs(m[r][i]):
            max_row = r
            max = abs(m[r][i])
    m[i], m[max_row] = m[max_row], m[i]


if __name__ == '__main__':
    # m = [[0,-2,6,-10], [-1,3,-6,5], [4,-12,8,12]]
    # m = [[1,-1,3,2], [3,-3,1,-1], [1,1,0,3]]
    m = [[0.5, 1.1, 3.1, 6], [5, 0.96, 6.5, 0.96], [2, 4.5, 0.36, 0.02]]
    print(gaussian_elimination_with_pivot(m))
    print(len(m))

    """  
    m = [[4,4,0,400], [-1,4,2,400], [0,-2,4,400]]   # aj Montri p80  [50, 50, 125]
    print(gaussian_elimination_with_pivot(m))
    """


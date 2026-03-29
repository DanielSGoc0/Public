import numpy as np
from scipy.special import comb

N = 20

T = [[0.0 for j in range(N)] for i in range(N)]

for i in range(N):
	T[i][0] = 1.0
	T[i][i] = 1.0
 
for i in range(1, N):
	for j in range(1, i):
		T[i][j] = np.sqrt(T[i - 1][j] + T[i - 1][j - 1])

T = np.array(T)
print(T)


for i in range(N):
	sum = 0.0
	for j in range(i + 1):
		# sum += (-1)**j * T[i][j] * comb(i, j)
		sum += (-1)**j * T[i][j]
	print(sum)

import numpy as np
import itertools as it

N = 5
D = 5

MIN = 10000000000.0
MIN_X = None
for attempt in range(10000000):
	X = np.random.normal(0.0, 1.0, (N, D))
	
	for k in range(N):
		X[k] /= np.sqrt(np.dot(X[k], X[k]))

	sum = 0.0
	for pi in it.permutations([i for i in range(N)]):
		cycles = 0
		visited = [False for k in range(N)]

		prod = 1.0
		for i in range(N):
			if visited[i]:
				continue
			cycles += 1

			x = i
			while not visited[x]:
				prod *= np.dot(X[x], X[pi[x]])
				visited[x] = True
				x = pi[x]

		sum += prod * 2.5**(N - cycles)

	if MIN > sum:
		MIN = sum
		MIN_X = X
	if attempt % 1000 == 0:
		print(attempt, "\t", sum, "\t", MIN)
		print(MIN_X)
		print(np.dot(MIN_X, MIN_X.T))

	# if sum < 1.0:
	# 	print(X)
	# 	break

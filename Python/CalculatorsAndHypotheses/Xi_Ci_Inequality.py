# Checks a nice little inequality.
# It works!
import numpy as np

N = 4

for k in range(1000000):

	C1 = np.random.random((N))
	C2 = np.random.random((N))
	X = np.random.random((N))

	WORKS = True

	for i in range(N - 1):
		# if X[i] / C1[i] < X[i + 1] / C1[i + 1]:
		# 	WORKS = False
		# 	break
		if X[i] / C2[i] < X[i + 1] / C2[i + 1]:
			WORKS = False
			break

	SUM1 = 0.0
	SUM2 = 0.0
	for i in range(N):
		SUM1 += C1[i]**2
		SUM2 += C2[i]**2

		if SUM1 > SUM2:
			WORKS = False
			break

	if WORKS and np.dot(X, C1) > np.dot(X, C2):
		print("SHIT")
		print(X)
		print(C1)
		print(C2)
		break

	if WORKS:
		print(k)

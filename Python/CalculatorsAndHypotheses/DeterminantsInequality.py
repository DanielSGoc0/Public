import numpy as np

N = 10

for k in range(1000):
	print(k)

	a = np.random.random((N + 1, N + 1))
	b = np.random.random((N + 1, N + 1))

	for j in range(N + 1):
		b[N][j] = 0.0

	A = a @ a.T
	B = b @ b.T

	CA = np.linalg.det(A) / np.linalg.det(A[:N, :N])
	CB = np.linalg.det(A + B) / np.linalg.det(A[:N, :N] + B[:N, :N])

	if CB < CA:
		print("SHIT")
		print(CB, CA)
		print(A)
		print(B)
		break

# Checks a hypothesis:
# If L is some integral wrt positive measure over Gaussian Process,
# And L' is its best approximation by finitely
# many evaluations, are the weights >= 0?
# NO
import numpy as np
from scipy.stats import ortho_group

for attempt in range(100000):
	if attempt % 100 == 0:
		print(attempt)

	N = 4
	n = 1

	D = np.diag(np.random.random((N)))
	U = ortho_group.rvs(dim=N)
	K = U @ D @ U.T
	m = np.random.random((N))
	# We consider L = m0 K0 + m1 K1 + ...

	BEST = 0.0
	BEST_w = None
	indices = [i for i in range(n)]
	while True:
		# print(K[indices, :][:, indices])
		a = K[indices, :] @ m
		w = np.linalg.solve(K[indices, :][:, indices], a)

		RES = a @ w
		if RES > BEST:
			BEST = RES
			BEST_w = w
		
		if indices[0] == N - n:
			break
		for i in range(n):
			if i == n - 1:
				indices[i] += 1
			elif indices[i] + 1 == indices[i + 1]:
				indices[i] = i
			else:
				indices[i] += 1
				break

	for k in range(n):
		if BEST_w[k] < -0.000001:
			print("DAMN")
			print(K)
			print(BEST)
			print(BEST_w)
			print(m)
			break



# Finds the optimal set of evaluations
# to approximate the integral of Brownian Bridge.
import numpy as np
import scipy

np.set_printoptions(suppress=True,
   formatter={'float_kind':'{:5.8f}\t'.format}, linewidth=200000)


def calc(dt_logs, p_logs, total_p):
	n = len(dt_logs)
	dt_logs[0] = 0.0
	for i in range(1, n):
		dt_logs[i] = dt_logs[1]
	# p_logs = np.zeros((n), dtype='double')

	# First reconstruct dt and p
	dt = np.zeros((n + 1), dtype='double')
	p = np.zeros((n), dtype='double')

	w = np.zeros((n + 1), dtype='double')
	w[n] = 1.0
	s = 1.0
	for i in range(n):
		w[i] = np.exp(dt_logs[i])
		s += w[i]
	for i in range(n + 1):
		dt[i] = w[i] / s
	
	s = 0.0
	for i in range(n):
		w[i] = np.exp(p_logs[i])
		s += w[i]
	for i in range(n):
		p[i] = total_p * w[i] / s

	# Prepare sequence t
	t = np.zeros((n + 2), dtype='double')
	t[0] = 0.0
	for i in range(0, n + 1):
		t[i + 1] = t[i] + dt[i]
	# so that t[i + 1] - t[i] = dt[i]

	# Next construct matrix Kinv
	K = np.zeros((n + 1, n + 1), dtype='double')

	for i in range(n):
		for j in range(n):
			K[i][j] = min(t[i + 1], t[j + 1]) - t[i + 1] * t[j + 1]
	for i in range(n):
		K[i][n] = t[i + 1] * (1.0 - t[i + 1]) / 2.0
		K[n][i] = t[i + 1] * (1.0 - t[i + 1]) / 2.0
	K[n][n] = 1.0/12.0

	# We can finally calculate the ANS
	ANS = K[n][n] - K[n, :n] @ np.linalg.solve(K[:n, :n] + np.diag(1.0 / p), K[:n, n])

	print()
	print(t)
	print(p)
	print(ANS)

	return ANS


def minimize():
	n = 4
	total_p = 1.0

	fun = lambda v: calc(v[0:n], v[n:(2*n)], total_p)
	# fun(np.array([np.random.normal(0.0, 1.0, (1)) for _ in range(2*n)]))
	start = np.array([np.random.normal(0.0, 1.0, (1)) for _ in range(2*n)])
	# for i in range(n, 2*n):
	# 	start[i] = 0.0

	RES = scipy.optimize.minimize(fun, start, options={"gtol": 0.0000000000000000000001})
	print("=================================================================")
	fun(RES.x)

minimize()

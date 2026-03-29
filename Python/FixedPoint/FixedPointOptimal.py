# Covariance is exp(-||X - Y||^2 / 2)
# This is an alternative, as simple as possible (and numerically stable)
# computation of results for fixed coordinates methods.
# Allows for error

import numpy as np
import scipy

np.set_printoptions(suppress=True,
   formatter={'float_kind':'{:5.8f}\t'.format}, linewidth=130)

err2 = 0.0
alpha = 0.8

def score(a):
	global err2, alpha
	# print(a)
	N = len(a)

	A = np.zeros((N + 1), dtype='double')
	B = np.zeros((N + 1), dtype='double')
	w = np.zeros((N + 1), dtype='double')

	# initial values
	w[0] = 1.0
	A[0] = 0.0 # half exp
	B[0] = 0.0 # no exp

	for n in range(N):
		w[n + 1] = 1 - np.exp(-a[n]**2) * (1 - err2 * w[n]/(err2 + w[n]))

		A[n + 1] = np.exp(-a[n]**2 / 2.0) * (A[n] + a[n] * w[n] / np.sqrt(err2 + w[n]))
		B[n + 1] = B[n] + a[n]**2

	print(B[N])
	print(A[N])
	print(alpha**2 - 2 * alpha * A[N] + B[N])
	return alpha**2 - 2 * alpha * A[N] + B[N]

N = 40
a0 = np.ones((N))
for i in range(N):
	# a0[i] = np.sqrt(8.0/(i + 8))
	a0[i] = 2.0/np.sqrt(i + 6)
	# a0[i] = 0.0001

# score(a0)

# score([0.5, 0.6, 0.7, 0.9, 0.3])
# score([0.5])

RES = scipy.optimize.minimize(score, a0, options={"gtol": 0.00000000001})
print("======================================================")
print(RES)
print(RES.x)



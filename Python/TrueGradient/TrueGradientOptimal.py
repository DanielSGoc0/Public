# Covariance is exp(-||X - Y||^2 / 2)
# This is an alternative, as simple as possible (and numerically stable)
# computation of results for fixed coordinates methods.
import numpy as np
import scipy

np.set_printoptions(suppress=True,
   formatter={'float_kind':'{:5.8f}\t'.format}, linewidth=130)

err2 = 0.0

def optimal(F0, w0, a0, N):
	a = np.zeros((N + 1), dtype='double')
	w = np.zeros((N + 1), dtype='double')
	F = np.zeros((N + 1), dtype='double')

	# initial values
	a[0] = a0
	w[0] = w0
	F[0] = F0

	for n in range(N):
		S = np.exp(-a[n]**2 / 2)
		w[n + 1] = np.sqrt(1 - S**2)
		F[n + 1] = S * (a[n] * w[n] + F[n])
		a[n + 1] = w[n + 1] / a[n] / S * (a[n] / S * F[n + 1] - w[n])

	return a


def score(a, F0, w0):
	global err2
	# print(a)
	N = len(a)

	w = np.zeros((N + 1), dtype='double')
	F = np.zeros((N + 1), dtype='double')

	# initial values
	w[0] = w0
	F[0] = F0

	for n in range(N):
		# calculate w[n + 1] and F[n + 1]
		w[n + 1] = 1 - np.exp(-a[n]**2) * (1 - err2 * w[n] / (err2 + w[n]))
		F[n + 1] = np.exp(-a[n]**2 / 2.0) * (F[n] + a[n] * w[n] / np.sqrt(err2 + w[n]))
		# w[n + 1] = 1 - np.exp(-a[n]**2) * (1 - w[n] / err2)
		# F[n + 1] = np.exp(-a[n]**2 / 2.0) * (F[n] + a[n] * w[n])

		# w[n + 1] = np.sqrt(1 - np.exp(-a[n]**2))
		# F[n + 1] = np.exp(-a[n]**2 / 2.0) * (F[n] + a[n] * w[n])

	print(F[N])

	return -F[N]

N = 30

F0 = 0.0
w0 = 1.0
# print(optimal(F0, 1.0, 1.0/F0, N))
a0 = np.ones((N + 1), dtype='double')
# a0 = optimal(F0, 1.0, 1.0/F0, N)

score(a0, F0, w0)

RES = scipy.optimize.minimize(lambda a: score(a, F0, w0), a0, options={"gtol": 0.00000000001})
print("======================================================")
print(RES)
print(RES.x)
print(a0)


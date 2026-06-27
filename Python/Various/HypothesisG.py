import numpy as np
import scipy

def H(r2, rho2):
	h2 = r2 / (1.0 - rho2 * np.exp(-r2))

	a2 = 0.0
	if h2 > 1.0:
		a2 = h2 + np.real(scipy.special.lambertw(-np.exp(-h2) * h2))

	res = np.sqrt(r2 * (1.0 - rho2 * np.exp(-r2)))
	res -= np.sqrt(a2 * (1.0 - np.exp(-a2)))
	res += 2.0 * (1.0 - np.exp(-a2 / 2.0))

	return res



for k in range(1):
	r2 = np.random.random() * 2.0
	rho2 = np.random.random()
	# x2 = np.random.random() * r2
	for i in range(100):
		x2 = (r2 * i) / 100.0

		LHS = np.exp(-r2 / 2.0) * np.sqrt(x2 * (1.0 - rho2)) + H(r2 - x2, np.exp(-x2))
		RHS = H(r2, rho2)

		print(RHS - LHS)

# Back when I tried to find the exact
# value of G function. But failed
import numpy as np
import scipy
import math

def G(r2, rho2):
	# h2 = r2 / (1.0 - rho2 * np.exp(-r2))
	# if h2 <= 1.0:
	# 	return np.sqrt(max(0.0, r2 * (1.0 - rho2 * np.exp(-r2))))
	# else:
	# 	r2bar = np.real(scipy.special.lambertw(max(0.000000001 - np.exp(-1.0), -np.exp(-h2) * h2))) + h2
	# 	if r2bar < 0.000001:
	# 		return r2 * (1.0 - r2bar/4.0 + r2bar**2/24.0)
	# 	return r2 * np.expm1(-r2bar/2.0)/(-r2bar/2.0)
	# w2 = np.real(-scipy.special.lambertw(-np.exp(-r2 - rho2 * np.exp(-r2)), 0) - rho2 * np.exp(-r2))
	# if 1.0 - np.exp(-r2) * rho2 > r2 or math.isnan(w2) or 1.0 - w2 < 0.0:
	# 	# return 0.0
	# 	return np.sqrt(r2 * (1.0 - np.exp(-r2) * rho2))
	# else:
	# 	# return 0.0
	# 	# return (1.0 + np.sqrt(1.0 - w)) * (1.0 - np.sqrt(rho2 * np.exp(-r2)))
	# 	return 2.0 * (1.0 - np.exp(-(r2 - w2)/2.0)) + w2
	# 	# return (1.0 + np.sqrt(1.0 - w)) * (1.0 - np.sqrt(rho2 * np.exp(-r2))) - np.sqrt(r2 * (1.0 - np.exp(-r2) * rho2))
	# 	# return (1.0 + np.sqrt(1.0 - w)) * (1.0 - np.sqrt(rho2 * np.exp(-r2))) - 2.0 * (1.0 - np.exp(-r2/2.0))


	a2 = rho2 * np.exp(-r2)
	h2 = r2/(1.0 - a2)
	w2 = np.real(scipy.special.lambertw(-h2 * np.exp(-h2), 0) + h2)
	if math.isnan(w2) or 1.0 - a2 > r2:
		return np.sqrt(r2 * (1.0 - np.exp(-r2) * rho2))
	else:
		return np.sqrt(r2 * (1.0 - np.exp(-r2) * rho2)) - np.sqrt(w2 * (1.0 - np.exp(-w2))) + 2.0 * (1.0 - np.exp(-w2/2.0))

r2 = 1.0
rho2 = 0.8

value = G(r2, rho2)
supremum = 0.0

N = 10000
for S in range(0, N + 1):
	x2 = (1.0 * S)/N * r2
	
	v = np.exp(-r2/2.0) * np.sqrt(x2 * (1.0 - rho2)) + G(r2 - x2, np.exp(-x2))
	supremum = max(supremum, v)

print(supremum, value)

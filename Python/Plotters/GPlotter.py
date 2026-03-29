# G functions plotter

import numpy as np
import matplotlib.pyplot as plt
import math
from matplotlib import cm
import scipy

K = 5
R2_MAX = 2.0

R2_AMT = 30
A_AMT = 100
R2 = [[(i * 1.0) / R2_AMT * R2_MAX for j in range(A_AMT + 1)] for i in range(R2_AMT + 1)]
A = [[(j * 1.0) / A_AMT for j in range(A_AMT + 1)] for i in range(R2_AMT + 1)]

X_AMT = 200
Y_AMT = 100
X = [[(i * 1.0) / X_AMT for j in range(Y_AMT + 1)] for i in range(X_AMT + 1)]
Y = [[(j * 1.0) / Y_AMT for j in range(Y_AMT + 1)] for i in range(X_AMT + 1)]

# G[0] is only zeros
# G[1] is filled below
G = [[[0.0 for j in range(A_AMT + 2)] for i in range(R2_AMT + 2)] for k in range(K + 1)]
G_x = [[[0.0 for j in range(A_AMT + 2)] for i in range(R2_AMT + 2)] for k in range(K + 1)]


def binsearch(f, y, a0, b0, N):
	a = a0
	b = b0
	for _ in range(N):
		h = (a + b)/2.0
		if f(h) < y:
			a = h
		else:
			b = h
	return (a + b)/2.0

if K >= 1:
	# for i  in range(R2_AMT + 1):
	# 	for j in range(A_AMT + 1):
	# 		G[1][i][j] = np.exp(-R2[i][j] / 2.0) * np.sqrt(R2[i][j] * A[i][j])
	for i  in range(R2_AMT + 1):
		for j in range(A_AMT + 1):
			if R2[i][j] + (1.0 - A[i][j]) * np.exp(-R2[i][j]) <= 1.0:
				G[1][i][j] = np.sqrt(R2[i][j] * (1.0 - np.exp(-R2[i][j]) * (1.0 - A[i][j])))
			else:
				G[1][i][j] = 2.0 * (1.0 - np.exp(-R2[i][j] / 2.0))

# if K >= 2:
# 	for i in range(R2_AMT + 1):
# 		for j in range(A_AMT + 1):
# 			z0 = 0.0
# 			if i > 0:
# 				func = lambda z: np.sqrt(4.0*z*np.exp(z) - 4.0*z)/(np.sqrt(A[i][j]) + np.sqrt(A[i][j] + 4.0*z*np.exp(z))) - np.sqrt(R2[i][j] - z)
# 				z0 = binsearch(func, 0.0, 0.0, R2[i][j], 50)
# 			v = np.exp(-R2[i][j] / 2.0) * (np.sqrt(z0 * A[i][j]) + np.sqrt(max(0.0, (R2[i][j] - z0) * (np.exp(z0) - 1))))
# 			G[2][i][j] = v

def best_G(r2, rho2):
	# return np.sqrt(r2 * (1.0 - np.exp(-r2) * rho2))
	# return 2.0/(1.0 + np.sqrt(rho2 * np.exp(-r2)))
	a2 = rho2 * np.exp(-r2)
	# w2 = np.real(-scipy.special.lambertw(-np.exp(-r2 - a2), 0) - a2)
	h2 = r2/(1.0 - a2)
	w2 = np.real(scipy.special.lambertw(-h2 * np.exp(-h2), 0) + h2)
	if math.isnan(w2) or 1.0 - a2 > r2:
		# return 0.0
		return np.sqrt(r2 * (1.0 - np.exp(-r2) * rho2))
		# return np.sqrt(2.0 * (1.0 - np.exp(-r2/2.0))*(1.0 - np.exp(-r2/2.0) + 1.0 - rho2 * np.exp(-r2/2.0))) - np.sqrt(r2 * (1.0 - np.exp(-r2) * rho2))
	else:
		# return w2/(1.0 - a2)
		# return (w2 + np.sqrt(2.0 * (r2 + a2 - 1.0)))**2/r2
		# return ((1.0 - a2 - w2)/(1.0 + np.sqrt(a2)))**2 + 1.0 - a2
		# return w2 + np.sqrt(2.0 * (r2 - w2)) - (1.0 - a2)
		# return w2 + np.sqrt(2.0 * (1.0 - w2 - a2)) - (1.0 - a2)
		# return (1.0 + np.sqrt(1.0 - w)) * (1.0 - np.sqrt(rho2 * np.exp(-r2)))
		# return (1.0 + np.sqrt(1.0 - w)) * (1.0 - np.sqrt(rho2 * np.exp(-r2))) - np.sqrt(r2 * (1.0 - np.exp(-r2) * rho2))
		# return (1.0 + np.sqrt(1.0 - w)) * (1.0 - np.sqrt(rho2 * np.exp(-r2))) - 2.0 * (1.0 - np.exp(-r2/2.0))
		# return np.sqrt(4.0 * (1.0 - np.exp(-r2/2.0))*(1.0 - np.sqrt(a2)))
		# return np.sqrt(2.0 * (1.0 - np.exp(-r2/2.0))*(1.0 - np.exp(-r2/2.0) + 1.0 - rho2 * np.exp(-r2/2.0))) - np.sqrt(r2 * (1.0 - np.exp(-r2) * rho2))
		# return np.sqrt(2.0 * (1.0 - np.exp(-r2/2.0))*(1.0 - np.exp(-r2/2.0) + 1.0 - rho2 * np.exp(-r2/2.0))) - np.sqrt(r2 * (1.0 - np.exp(-r2) * rho2))
		# return np.sqrt(4.0 * (1.0 - np.exp(-r2/2.0))/(1.0 + np.sqrt(a2))) * (1.0 - a2)
		# return 1 - a2
		# return 2.0 * np.sqrt(2.0 / (1.0/(1.0 - np.exp(-r2/2.0))**2 + (1.0 - np.sqrt(a2))/(1.0 - np.sqrt(a2))**2))
		# return 2.0 * (1.0 - np.sqrt(a2))

		return np.sqrt(r2 * (1.0 - np.exp(-r2) * rho2)) - np.sqrt(w2 * (1.0 - np.exp(-w2))) + 2.0 * (1.0 - np.exp(-w2/2.0))


		
def best_G2(r2, v2):
	w2 = r2 * v2
	rho2 = np.exp(w2) - w2 * np.exp(r2)
	return best_G(r2, rho2)
	

def approx_G(R2, A, k):
	global G, R2_MAX, R2_AMT, A_AMT

	R2 = R2 / R2_MAX * R2_AMT
	A = A * A_AMT

	i = math.floor(R2)
	ci = R2 - i
	j = math.floor(A)
	cj = A - j

	v = 0.0
	v += G[k][i][j] * (1.0 - ci) * (1.0 - cj)
	v += G[k][i + 1][j] * ci * (1.0 - cj)
	v += G[k][i][j + 1] * (1.0 - ci) * cj
	v += G[k][i + 1][j + 1] * ci * cj

	return v

def approx_G2(r2, v2, k):
	w2 = r2 * v2
	rho2 = np.exp(w2) - w2 * np.exp(r2)
	return approx_G(r2, rho2, k)

def approx_Gx(R2, A, k):
	global G_x, R2_MAX, R2_AMT, A_AMT

	R2 = R2 / R2_MAX * R2_AMT
	A = A * A_AMT

	i = math.floor(R2)
	ci = R2 - i
	j = math.floor(A)
	cj = A - j

	v = 0.0
	v += G_x[k][i][j] * (1.0 - ci) * (1.0 - cj)
	v += G_x[k][i + 1][j] * ci * (1.0 - cj)
	v += G_x[k][i][j + 1] * (1.0 - ci) * cj
	v += G_x[k][i + 1][j + 1] * ci * cj

	return v

def f(R2, A, x, y, k):
	global G, R2_MAX, R2_AMT, A_AMT

	x *= R2
	y *= np.sqrt(A)
	
	R2_ = max(0.0, R2 - x)
	A_ = max(0.0, 1.0 - np.exp(-x) * (1.0 - A + y**2))

	v = approx_G(R2_, A_, k)
	v += np.exp(-R2 / 2.0)*np.sqrt(x)*y

	return v

if K >= 2:
	for k in range(2, K + 1):
		for i in range(R2_AMT + 1):
			for j in range(A_AMT + 1):
				G_x[k][i][j] = 0.0

				if R2[i][j] + (1.0 - A[i][j]) * np.exp(-R2[i][j]) <= 1.0:
					G[k][i][j] = np.sqrt(R2[i][j] * (1.0 - np.exp(-R2[i][j]) * (1.0 - A[i][j])))
					continue

				supremum = f(R2[i][j], A[i][j], X[X_AMT][Y_AMT], Y[X_AMT][Y_AMT], k - 1)
				supremum_at = (X[X_AMT][Y_AMT], Y[X_AMT][Y_AMT])

				for a in range(X_AMT + 1):
					v = f(R2[i][j], A[i][j], X[a][Y_AMT], 1.0, k - 1)
					if v > supremum:
						supremum = v
						supremum_at = R2[i][j] * X[a][Y_AMT]

				G[k][i][j] = supremum

				r2 = R2[i][j] - supremum_at
				rho2 = np.exp(-supremum_at)
				# G_x[k][i][j] = r2 / (1.0 - np.exp(-R2[i][j]))

				a2 = rho2 * np.exp(-r2)
				w2 = np.real(-scipy.special.lambertw(-np.exp(-r2 - a2), 0) - a2)
				if math.isnan(w2):
					G_x[k][i][j] = 0.0
				else:
					G_x[k][i][j] = w2/r2
				# if r2 + rho2 * np.exp(-r2) <= 1.0:
				# 	G_x[k][i][j] = r2
				# 	# G_y[k][i][j] = rho2
				# 	G_y[k][i][j] = 1.0

# ========================================================================================================================

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
ax.autoscale(enable=None, axis="x", tight=True)
ax.set_xlabel('r2')
ax.set_ylabel('rho2')
ax.set_zlabel('value')

# Z = np.array([[G[K][i][j] for j in range(A_AMT + 1)] for i in range(R2_AMT + 1)])

Z = np.array([[0.0 for j in range(A_AMT + 1)] for i in range(R2_AMT + 1)])
# We set T = X[i][j]  and   t = X[i][j] * Y[i][j]
for i in range(1, R2_AMT + 1):
	for j in range(A_AMT + 1):
		# r2 = R2[i][j]
		# # # a2 = R2[i][j] / R2_MAX
		# v2 = A[i][j]
		# w2 = r2 * A[i][j]
		# # # r2 = w2 - np.log(w2 + a2)
		# rho2 = np.exp(w2) - w2 * np.exp(r2)
		# if rho2 < 0.0 or rho2 > 1.0 or w2 > r2 or r2 < 0.0 or r2 > R2_MAX or rho2 * np.exp(-r2) + r2 < 1.0:
		# 	Z[i][j] = 0.0
		# 	continue
		# # Z[i][j] = approx_G(r2, 1.0 - rho2, K) - (1.0 - v2) * 2.0 * (1.0 - np.exp(-r2/2.0)) - v2 * r2
		# # Z[i][j] = approx_G2(r2, v2, K)
		# # Z[i][j] = approx_G2(r2, v2, K) - 2.0 * (1.0 - np.exp(-r2/2.0))
		# # Z[i][j] = approx_G(r2, 1.0 - rho2, K)
		# # Z[i][j] = approx_Gx(r2, 1.0 - rho2, K)
		# # Z[i][j] = best_G(r2, rho2)



		# S2 = np.exp(1.0) * (A[i][j] + 1.0)
		# T = R2[i][j] / R2_MAX
		# h2 = -scipy.special.lambertw(min(max(-np.exp(-1.0) + 0.00000001, -1.0/S2), -0.00000001), -1)
		# h2 = np.real(h2)
		# r2 = h2 - T
		# rho2 = S2 * np.exp(-T) * T
		# if r2 < 0.0 or rho2 > 1.0:
		# 	Z[i][j] = 0.0
		# 	continue

		r2 = R2[i][j]
		rho2 = A[i][j]
		# Z[i][j] = best_G(r2, rho2)
		Z[i][j] = best_G(r2, rho2) - approx_G(r2, 1.0 - rho2, K)
		# h2 = r2 / (1.0 - rho2 * np.exp(-r2))

		# if h2 < 1.0:
		# 	Z[i][j] = 1.0
		# 	continue

		# r2 = R2[i][j]
		# h2 = 1.0 + A[i][j] * (R2_MAX - 1.0)
		# rho2 = np.exp(r2) * (1.0 - r2 / h2)

		# r2_of_h2 = scipy.special.lambertw(max(0.00000001 - np.exp(-1.0), -np.exp(-h2) * h2)) + h2

		# if rho2 > 1.0 or rho2 < 0.0:
		# 	Z[i][j] = 0.0
		# 	continue

		# if Y[i][j] > X[i][j] or 1.0 - np.exp(T-t) * t/T < 0.0:
		# 	Z[i][j] = 0.0
		# else:
		# 	# Z[i][j] = np.exp(-t/2.0) * approx_G(T-t, np.exp(T-t) * t/T, K) + np.exp((T - t)/2.0) * (t + 1)/np.sqrt(T)
		# 	# Z[i][j] = approx_G(T-t, 1.0 - np.exp(T-t) * t/T, K)/np.sqrt(T) - (T - t)/T
		# 	# Z[i][j] = best_G(T-t, np.exp(T-t) * t/T)/np.sqrt(T) - (T - t)/T
		# 	# Z[i][j] = best_G(X[i][j], Y[i][j])
		# Z[i][j] = np.log(approx_G(r2, 1.0 - rho2, K)/np.sqrt(r2))
		# Z[i][j] = np.log(approx_G(r2, 1.0 - rho2, K)**2/r2 + np.exp(-r2) * rho2)
		# Z[i][j] = best_G(r2, rho2)
		# Z[i][j] = approx_G(r2, 1.0 - rho2, K) - np.sqrt(r2 * (1.0 - np.exp(-r2) * rho2))
		# Z[i][j] = approx_G(r2, 1.0 - rho2, K) - approx_G(r2, 1.0 - rho2, 1)
		# Z[i][j] = approx_G(r2, 1.0 - rho2, K) / (1.0 - np.sqrt(rho2) * np.exp(-r2/2.0)) / (approx_G(1.0, 1.0 - rho2, K) / (1.0 - np.sqrt(rho2) * np.exp(-1.0/2.0)))
		# Z[i][j] = approx_G(r2, 1.0 - rho2, K) / (r2 / h2) * (1.0 + np.exp(-r2_of_h2/2.0))/2.0
		# Z[i][j] = approx_Gx(r2, 1.0 - rho2, K)
		# Z[i][j] = (best_G(r2, rho2)**2 - approx_G(r2, 1.0 - rho2, K)**2)/r2
		# if r2 <= 1.0:
		# 	Z[i][j] = 0.0
		# if r2 + rho2 * np.exp(-r2) <= 1.0:
		# 	print("wtf???", 1.0 - (r2 + rho2 * np.exp(-r2)))
		# 	Z[i][j] = -10.0

		# Z[i][j] = (approx_G(r2, 1.0 - rho2, K)**2 - approx_G(h2, 1.0 - 0.0, K)**2)/T - (approx_G(0, 1.0 - 1.0, K)**2 - approx_G(h2, 1.0 - 0.0, K)**2)/h2
		# Z[i][j] = best_G(r2, rho2) - approx_G(r2, 1.0 - rho2, K)
		# Z[i][j] = best_G(r2, rho2)




# Z2 = np.array([[0.0 for j in range(A_AMT + 1)] for i in range(R2_AMT + 1)])
# for i in range(1, R2_AMT + 1):
# 	for j in range(0, A_AMT + 1):
# 		if j == 0 or j == A_AMT:
# 			Z2[i][j] = 0.0
# 		else:
# 			Z2[i][j] = np.log((Z[i][j - 1] + Z[i][j + 1] - 2.0 * Z[i][j]) * A_AMT**2)

surf = ax.plot_surface(R2, A, Z, cmap=cm.coolwarm, linewidth=0, antialiased=False)
# surf = ax.plot_surface(R2, A, Z2, cmap=cm.coolwarm, linewidth=0, antialiased=False)
# surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm, linewidth=0, antialiased=False)

plt.show()

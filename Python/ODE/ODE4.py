import numpy as np
import scipy
import math

S2 = math.sqrt(2)
LS2 = math.log(2)/2.0

def get_v(N, X):
	v = []
	for i in range(N):
		for j in range(i):
			v.append(X[i][j])
	return v


def get_X(N, v):
	X = np.zeros((N, N), dtype='double')
	k = 0
	for i in range(N):
		for j in range(i):
			X[i][j] = v[k]
			k += 1
	return X


def make_vector(funs, x, k):
	res = np.zeros((k), dtype='double')
	for i in range(k):
		res[i] = funs[i](x)
	return res


# ===============================================================================

def create_callable_integral(f, a, b, nodes):
	xs = []
	ys = []
	dys = []

	xs.append(0)
	ys.append(0)
	dys.append(f(0))
	for i in range(1, nodes + 1):
		x = (a + b)/2.0 - (b - a)/2.0 * math.cos(math.pi * i / nodes)
		xs.append(x)
		ys.append(ys[-1] + scipy.integrate.fixed_quad(lambda v: [f(x) for x in v], xs[-2], xs[-1], n=5)[0])
		dys.append(f(x))

	return scipy.interpolate.CubicHermiteSpline(xs, ys, dys)


def calc(X, R):
	global S2, LS2

	X = np.array(X, dtype='double')
	N = X.shape[0]

	M = np.zeros((N, N), dtype='double')
	for i in range(N):
		for j in range(N):
			M[i][j] = math.exp(np.dot(X[i], X[j]))
	V = np.linalg.cholesky(M)

	dANS = np.zeros((N, N), dtype='double')
	ANS = -2*math.log(1 - R*R)

	# find sequence of functions
	gs = [lambda r: S2 * r]
	
	for k in range(1, N):
		a = scipy.linalg.solve_triangular(np.transpose(V), X[k], lower=False)
		gk = lambda r: math.exp(LS2 + sum([a[i] * gs[i](r) for i in range(k)]))
		gs.append(create_callable_integral(gk, 0, R, 100))

	g = make_vector(gs, R, N)
	x = scipy.linalg.solve_triangular(V, g, lower=True)
	ANS -= np.dot(x, x)

	# calculate all the derivatives
	# for i in range(N):
	# 	for j in range(i):
	# 		# dV is the derivative of V inverse transpose
	# 		dV = np.zeros((N, N), dtype='double')
	# 		dgs = [lambda r: 0]

	# 		# calculate dV
	# 		v1 = np.zeros(N, dtype='double')
	# 		v1[i] = 1
	# 		v1 = scipy.linalg.solve_triangular(V, v1, lower=True)

	# 		v2 = np.zeros(N, dtype='double')
	# 		for k in range(N):
	# 			v2[k] = M[k][i] * X[k][j]
	# 		v2 = scipy.linalg.solve_triangular(V, v2, lower=True)

	# 		for k in range(N):
	# 			for l in range(k):
	# 				dV[k][l] = v1[k] * v2[l] + v1[l] * v2[k]
	# 			dV[k][k] = v1[k] * v2[k]
	# 		dV = scipy.linalg.solve_triangular(np.transpose(V), -np.transpose(dV), lower=False)
		
	# 		# calculate sequence dgs. With lower precision
	# 		for k in range(1, N):
	# 			da = np.dot(dV, X[k])
	# 			if i == k:
	# 				v3 = np.zeros(N, dtype='double')
	# 				v3[j] = 1
	# 				v3 = scipy.linalg.solve_triangular(np.transpose(V), v3, lower=False)
	# 				da += v3
	# 			a = scipy.linalg.solve_triangular(np.transpose(V), X[k], lower=False)

	# 			dgk = lambda r: S2 * (sum([da[i] * gs[i](r) + a[i] * dgs[i](r) for i in range(k)])) * math.exp(sum([a[i] * gs[i](r) for i in range(k)]))
	# 			dgs.append(create_callable_integral(dgk, 0, R, 20))

	# 		# calculate the final derivative
	# 		dg = make_vector(dgs, R, N)
	# 		dx = scipy.linalg.solve_triangular(V, dg, lower=True) + np.dot(np.transpose(dV), g)

	# 		dANS[i][j] = -2*np.dot(x, dx)

	print(ANS)
	# print(dANS)
	# return (ANS, get_v(N, dANS))

	return ANS


# =========================================================================================

def initial_guess(N, R):
	global S2

	X = np.zeros((N, N), dtype='double')
	funs = []
	funs.append(lambda r: S2 * r)

	for k in range(1, N):

		M = np.zeros((k, k), dtype='double')
		for i in range(k):
			for j in range(k):
				M[i][j] = math.exp(np.dot(X[i], X[j]))
		V = np.linalg.cholesky(M)

		g = make_vector(funs, 1.0 * R * k/N, k)
		b = scipy.linalg.solve_triangular(V, g, lower=True)
		for i in range(k):
			X[k][i] = b[i]

		if k == N - 1:
			break

		a = scipy.linalg.solve_triangular(np.transpose(V), b, lower=False)
		f = lambda r: S2 * math.exp(sum([a[i] * funs[i](r) for i in range(k)]))
		funs.append(create_callable_integral(f, 0, R, 100))

	return X




def minimize(N, R):
	v0 = initial_guess(N, R)
	# print(v0)

	fun = lambda v: calc(get_X(N, v), R)

	# print(scipy.optimize.minimize(fun, get_v(N, v0), jac=True))
	RES = scipy.optimize.minimize(fun, get_v(N, v0), jac=False)
	X = get_X(N, RES.x)
	print(RES)
	for i in range(N - 1):
		print(np.dot(X[i + 1] - X[i], X[i + 1] - X[i]))



# v = [-1.0, 1.0, 0.01]
# print(error(get_X(3, v), 0.1))
# print(error(get_X(3, v), 0.2))
# print(error(get_X(3, v), 0.3))
# print(error(get_X(3, v), 0.4))

# print(error([[0.0, 0.0], [1.0, 0.0]], 1))
# print(error(get_X(5, [0.2641, 0.5078, 0.1668, 0.7518, 0.3062, 0.1462, 0.9371, 0.4536, 0.2039, 0.1123]), 1))
# print(error(get_X(3, [0.5229, 0.7670, 0.3853]), 1))

minimize(10, 0.9)

# h = 0.001
# (ANS1, dANS1) = calc([[0.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0], [0.5, 0.3, 0.0, 0.0], [0.1, 0.3, 0.5, 0.0]], 0.8)
# (ANS2, dANS2) = calc([[0.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0], [0.5, 0.3, 0.0, 0.0], [0.1, 0.3, 0.5, 0.0]], 0.8)

# print((ANS2 - ANS1)/h)
# print(dANS1)

# calc([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 0.3, 0.0]], 0.8)

# print(calc([[0.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0], [0.5, 0.3, 0.0, 0.0], [0.1, 0.3, 0.5, 0.0]], 0.8))

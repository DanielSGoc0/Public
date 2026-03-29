# This program is useful if we don't search for optimal nodes with small error (like 10^(-4))
# Instead, it is useful if the error is relatively big, say >= 1.
# Added the optional measurement error.
# changed the covariance function from e^(-x^2/2) to e^(-x^2), which makes things MUCH cleaner.

import numpy as np
import scipy
import math

err2 = 0.00

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
		x = a + (b - a) * i / nodes
		xs.append(x)
		ys.append(ys[-1] + scipy.integrate.fixed_quad(lambda v: [f(x) for x in v], xs[-2], xs[-1], n=5)[0])
		dys.append(f(x))

	return scipy.interpolate.CubicHermiteSpline(xs, ys, dys)


def calc(X, T, h):
	global err2
	P = 2*T + math.log((1 + math.exp(-2*T))/2.0)
	STEPS = math.ceil(P/h)

	X = np.array(X, dtype='double')
	N = X.shape[0]

	M = np.zeros((N, N), dtype='double')
	for i in range(N):
		for j in range(N):
			M[i][j] = math.exp(-np.dot(X[i] - X[j], X[i] - X[j]))
		M[i][i] += err2
	print(X)
	V = np.linalg.cholesky(M)

	ANS = P - math.log(2.0 - math.exp(-P))

	# find sequence of functions
	gs = [lambda p: 1 - math.exp(-p)]
	L2 = [np.dot(X[i], X[i]) for i in range(N)]

	for k in range(1, N):
		a = scipy.linalg.solve_triangular(np.transpose(V), X[k], lower=False)
		print(a)
		print(L2[k])
		gk = lambda p: math.exp(-p + 2*sum([a[i] * gs[i](p) for i in range(k)]) - L2[k])
		gs.append(create_callable_integral(gk, 0.0, P, STEPS))

	g = make_vector(gs, P, N)
	x = scipy.linalg.solve_triangular(V, g, lower=True)
	ANS -= np.dot(x, x)
	print("opt: ", x)

	print(ANS)

	print(X)

	return ANS


# =========================================================================================

def initial_guess(N, T, h):
	global err2
	P = 2*T + math.log((1 + math.exp(-2*T))/2.0)
	STEPS = math.ceil(P/h)

	X = np.zeros((N, N), dtype='double')
	gs = [lambda p: 1 - math.exp(-p)]

	for k in range(1, N):

		M = np.zeros((k, k), dtype='double')
		for i in range(k):
			for j in range(k):
				M[i][j] = math.exp(-np.dot(X[i] - X[j], X[i] - X[j]))
			M[i][i] += err2
		V = np.linalg.cholesky(M)

		g = make_vector(gs, 1.0 * P * k/N, k)
		b = scipy.linalg.solve_triangular(V, g, lower=True)
		for i in range(k):
			X[k][i] = b[i]

		if k == N - 1:
			break

		L2 = np.dot(b, b)
		a = scipy.linalg.solve_triangular(np.transpose(V), b, lower=False)
		gk = lambda p: math.exp(-p + 2*sum([a[i] * gs[i](p) for i in range(k)]) - L2)
		gs.append(create_callable_integral(gk, 0.0, P, STEPS))

	return X



def minimize(N, T, h):
	v0 = initial_guess(N, T, h)
	print(v0)

	fun = lambda v: calc(get_X(N, v), T, h)

	RES = scipy.optimize.minimize(fun, get_v(N, v0))
	X = get_X(N, RES.x)
	print(RES)
	for i in range(N - 1):
		print(np.dot(X[i + 1] - X[i], X[i + 1] - X[i]))


# print(calc([[0.0, 0.0, 0.0], [0.4385, 0.0, 0.0], [0.6541, 0.2786, 0.0]], 10, 0.01))

minimize(6, 30, 0.01)

# v0 = [0.5895, 0.8148, 0.6386]
# calc(get_X(3, v0), 30, 0.01)



# err2 = 0, T = 30

# N = 2
# 56.77358213537121
#  [[0.        	0.        ]
# 	[0.71784966	0.        ]]

# N = 3
# 56.010998339281045
#  [[0.        	0.        	0.        ]
# 	[0.59408316	0.        	0.        ]
# 	[0.82038922	0.64787573	0.        ]]

# N = 4
# 55.30103866894774
#  [[0.        	0.        	0.        	0.        ]
# 	[0.51391845	0.        	0.        	0.        ]
# 	[0.72843289	0.51999086	0.        	0.        ]
# 	[0.9019932  0.72391884	0.62842214	0.        ]]

# N = 5
# 54.63288391830544
#  [[0.         0.         	0.         	0.         	0.        ]
#	[0.45281409 0.         	0.         	0.         	0.        ]
#	[0.66513222 0.43504628 	0.         	0.         	0.        ]
#	[0.84431199 0.62413391 	0.50299547 	0.         	0.        ]
#	[0.95244294 0.78392291 	0.6929246  	0.62226584 	0.        ]]

# N = 6
# 54.000246946419416
#  [[0.         0.         	0.         	0.         	0.         	0.        ]
#	[0.40348383 0.         	0.         	0.         	0.         	0.        ]
#	[0.61532657 0.37018834 	0.         	0.         	0.         	0.        ]
#	[0.80108856 0.55283418 	0.41802513 	0.         	0.         	0.        ]
#	[0.92175496 0.7141362  	0.5905098  	0.50076833 	0.         	0.        ]
#	[0.97648082 0.81875835 	0.74374036 	0.68759735 	0.61736459 	0.        ]]

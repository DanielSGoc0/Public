# kernel is: <f(x0), f(x1)> = exp(-||x0 - x1||^2)
# Let M(x) be the solution of ODE y(0) = x and y'(t) = f(y(t)) and M(x) = y(T)
# Then we wish to calculate M^{-1}(0) using evaluations of M(x) at adaptovely chosen points.

import scipy
import numpy as np
import random
from pprint import pp

import scipy.integrate


N = 6
T = 1.0
# A is strictly lower diagonal
A = [[0.0 for j in range(N + 1)] for i in range(N + 1)]

# ==============================     SETTINGS AND HELPERS    ==================================

EPS = 1e-200
INT_ORDER = 7
SPLINE_ORDER = 19
OFFSET = (N + 2)*(SPLINE_ORDER + 2)
# OFFSET = 100
h = T / 2**INT_ORDER
original_range = range(-OFFSET, 2**INT_ORDER + 1 + OFFSET)
nodes = [o * h for o in original_range]
processed_nodes = [nodes[0] for _ in range(SPLINE_ORDER + 1)] + nodes[(SPLINE_ORDER+1)//2:-(SPLINE_ORDER+1)//2] + [nodes[-1] for _ in range(SPLINE_ORDER + 1)]


def logcosh(t):
	t = abs(t)
	return t + np.log((1 + np.exp(-2*t))/2.0)

def neg(y):
	for i in range(len(y)):
		y[i] = -y[i]
	return y

# Returns the Bspline corresponding to f(t)
def make_BSpline(f):
	global EPS, INT_ORDER, SPLINE_ORDER, OFFSET, h, original_range, nodes

	vals = [f(n) for n in nodes]
	return scipy.interpolate.make_interp_spline(nodes, vals, k=SPLINE_ORDER)

# Returns the Bspline corresponding to the array y
def make_BSpline_tab(y):
	global EPS, INT_ORDER, SPLINE_ORDER, OFFSET, h, original_range, nodes
	return scipy.interpolate.make_interp_spline(nodes, y, k=SPLINE_ORDER)

# Sum up BSplines with weights, with indexes in range(m).
def sum_BSplines(weights, Bsplines, m):
	global EPS, INT_ORDER, SPLINE_ORDER, OFFSET, h, original_range, nodes, processed_nodes

	c = [0.0 for o in original_range]
	for o in original_range:
		for k in range(m):
			c[o + OFFSET] += weights[k] * Bsplines[k].c[o + OFFSET]
	
	return scipy.interpolate.BSpline(processed_nodes, c, k=SPLINE_ORDER)

# solve an ODE y'(t) = f(t, y(t)), y(0) = y0 for t in nodes.
def solve_ODE(f, y0):
	global EPS, INT_ORDER, SPLINE_ORDER, OFFSET, h, original_range, nodes

	M = len(y0)

	sol1 = scipy.integrate.solve_ivp(f, (0, T + OFFSET * h), y0, method='DOP853', t_eval=nodes[OFFSET:], first_step=h, max_step=h)
	sol2 = scipy.integrate.solve_ivp(lambda t, z: neg(f(-t, z)), (0, OFFSET * h), y0, method='DOP853', t_eval=nodes[OFFSET + 1:2*OFFSET + 1], first_step=h, max_step=h)

	RES = []

	for m in range(M):
		y = []
		for i in original_range:
			if i < 0:
				y.append(sol2.y[m][-i - 1])
			else:
				y.append(sol1.y[m][i])
		RES.append(y)
		
	return RES


# ==================================================================================


# L[i][j](ti, tj) = <yi(ti), yj(tj)>
L_0_ = [[None for j in range(N + 1)] for i in range(N + 1)]
L_T_ = [[None for j in range(N + 1)] for i in range(N + 1)]
L__0 = [[None for j in range(N + 1)] for i in range(N + 1)]
L__T = [[None for j in range(N + 1)] for i in range(N + 1)]

L_00 = [[None for j in range(N + 1)] for i in range(N + 1)]
L_T0 = [[None for j in range(N + 1)] for i in range(N + 1)]
L_0T = [[None for j in range(N + 1)] for i in range(N + 1)]
L_TT = [[None for j in range(N + 1)] for i in range(N + 1)]

# L2[i](ti) = <yi(ti), yi(ti)>
L2__ = [None for i in range(N + 1)]
L2_0 = [None for i in range(N + 1)]
L2_T = [None for i in range(N + 1)]


def get_V_0T(i, j):
	global EPS, INT_ORDER, SPLINE_ORDER, OFFSET, h, original_range, nodes, L2__, L_0_

	y = [np.exp(2*L_0_[i][j](n) - L2__[j](n)) for n in nodes[OFFSET:-OFFSET]]
	return np.log(max(scipy.integrate.romb(y, h), EPS))


# (dV[i][0]/dti, ... , dV[i][i - 1]/dti, dL2[i]/dti)
def get_V_derivative(t, V):
	global A

	i = len(V) - 1
	
	# S[i][j](ti, tj) = <f(yi(ti)), yj(tj)>
	S__0 = [0.0 for j in range(i + 1)]
	S__T = [0.0 for j in range(i)]
	
	for j in range(0, i + 1):
		for k in range(j):
			S__0[j] += A[j][k] * S__T[k]
		if j != i:
			S__T[j] = S__0[j] + np.exp(-V[-1] + V[j])
	
	RES = [S__0[j] + S__T[j] for j in range(i)]
	RES.append(2*np.tanh(t) + 2*S__0[i])

	return RES


def extract_L(ODE_solutions):
	global EPS, INT_ORDER, SPLINE_ORDER, OFFSET, h, original_range, nodes, L__0, L__T, L_0_, L_T_, L_00, L_0T, L_T0, L_TT

	i = len(ODE_solutions) - 1

	S__0_integrated = [[0.0 for o in original_range] for j in range(i)]
	S__T_integrated = [[0.0 for o in original_range] for j in range(i)]

	for o in original_range:
		for j in range(i):
			for k in range(j):
				S__0_integrated[j][o + OFFSET] += A[j][k] * S__T_integrated[k][o + OFFSET]
			S__T_integrated[j][o + OFFSET] = ODE_solutions[j][o + OFFSET] - S__0_integrated[j][o + OFFSET]

	y0 = [[S__0_integrated[j][o + OFFSET] - S__0_integrated[j][OFFSET] + L_00[i][j] for o in original_range] for j in range(i)]
	yT = [[S__T_integrated[j][o + OFFSET] - S__T_integrated[j][OFFSET] + L_0T[i][j] for o in original_range] for j in range(i)]

	for j in range(i):
		L__0[i][j] = make_BSpline_tab(y0[j])
		L__T[i][j] = make_BSpline_tab(yT[j])
		L_0_[j][i] = L__0[i][j]
		L_T_[j][i] = L__T[i][j]

		L_T0[i][j] = L__0[i][j](T)
		L_TT[i][j] = L__T[i][j](T)
		L_0T[j][i] = L_T0[i][j]
		L_TT[j][i] = L_TT[i][j]



def make_L_T_(i, j):
	global EPS, INT_ORDER, SPLINE_ORDER, OFFSET, h, original_range, nodes, L_00, L2__, L__0, L_0T

	# We consider y(tj) = V[j][i](tj, T) - 2 L[j][i](tj, 0)
	y0 = get_V_0T(j, i) - 2*L_00[j][i]
	ODE_solution = solve_ODE(lambda t, y: [np.exp(-L2__[j](t) + y[0] + 2*L__0[j][i](t))], [y0])[0]

	for o in original_range:
		n = nodes[o + OFFSET]
		ODE_solution[o + OFFSET] += L__0[j][i](n) - y0 + L_0T[j][i] - L_00[j][i]
	
	L_T_[i][j] = make_BSpline_tab(ODE_solution)
	L__T[j][i] = L_T_[i][j]


def calc():
	global N, T, A

	# print(A)

	# At i-th stage, we assume the knowledge of all the functions and variables for i', j' < i
	for i in range(0, N + 1):
		# Derive L_0_[i][j] for j < i
		for j in range(i):
			# <yi(0), yj(tj)> = sum_{k = 0}^{i - 1} A_{i, k} * <yj(tj), yk(T)>
			L_0_[i][j] = sum_BSplines(A[i], L__T[j], i)
			L__0[j][i] = L_0_[i][j]

			L_00[i][j] = L_0_[i][j](0)
			L_00[j][i] = L_00[i][j]
			L_0T[i][j] = L_0_[i][j](T)
			L_T0[j][i] = L_0T[i][j]

		# Derive L2_0[i]
		L2_0[i] = 0.0
		for k in range(i):
			for l in range(i):
				L2_0[i] += A[i][k] * A[i][l] * L_TT[k][l]
		L_00[i][i] = L2_0[i]

		# Next, derive L__0[i][j] and L__T[i][j] for j <= i, as well as L2__[i]
		# Start by finding V[i][j](0, T) = ln(int_{0}^{T} exp(2*L[i][j](0, t) - L2[j](t)) dt)
		y0 = [get_V_0T(i, j) for j in range(i)]
		y0.append(L2_0[i])

		ODE_solutions = solve_ODE(get_V_derivative, y0)

		# From ODE_solutions we extract L2[i](ti)
		L2__[i] = make_BSpline_tab(ODE_solutions[i])
		L2_T[i] = L2__[i](T)

		# Next from ODE_solutions we extract L__0[i][j] and L__T[i][j] for j < i.
		extract_L(ODE_solutions)

		# And then, extract L[i][i].
		# L = (L2__[i](ti) + L2__[i](tj)) / 2.0 - logcosh(ti - tj)
		L_0_[i][i] = make_BSpline(lambda tj: (L2_0[i] + L2__[i](tj)) / 2.0 - logcosh(tj))
		L_T_[i][i] = make_BSpline(lambda tj: (L2_T[i] + L2__[i](tj)) / 2.0 - logcosh(T - tj))
		L__0[i][i] = L_0_[i][i]
		L__T[i][i] = L_T_[i][i]

		L_T0[i][i] = (L2_0[i] + L2_T[i]) / 2.0 - logcosh(T)
		L_0T[i][i] = L_T0[i][i]
		L_TT[i][i] = L2_T[i]

		# It remains to find L_T_[i][j] for all j < i.
		for j in range(i):
			make_L_T_(i, j)
			# This is meant to give some idea about the error introduced at each step.
			# print("error: (", i, " ", j, ") = ", L_TT[i][j] - L_T_[i][j](T))

		# And update the fixed values too.
		for j in range(i):
			L_00[j][i] = L_00[i][j]
			L_0T[j][i] = L_T0[i][j]
			L_T0[j][i] = L_0T[i][j]
			L_TT[j][i] = L_TT[i][j]

	for i in range(N + 1):
		for j in range(N + 1):
			print('%.4f'%float(L_00[i][j]), end='\t')
			print('%.4f'%float(L_0T[i][j]), end='\t')
		print()
		for j in range(N + 1):
			print('%.4f'%float(L_T0[i][j]), end='\t')
			print('%.4f'%float(L_TT[i][j]), end='\t')
		print()
	# Return ||yN(T)||^2
	# print(L2_T[N])
	return L2_T[N]


# ==================================================================================

def get_v():
	global N, A

	v = []
	for i in range(N + 1):
		for j in range(i):
			v.append(A[i][j])
	return v

def set_A(v):
	global A, N

	A = np.zeros((N + 1, N + 1), dtype='double')
	k = 0
	for i in range(N + 1):
		for j in range(i):
			A[i][j] = v[k]
			k += 1

def set_default_A():
	global A, N
	for i in range(N + 1):
		for j in range(N + 1):
			if j < i:
				A[i][j] = -1.0
			else:
				A[i][j] = 0.0

def set_random_A():
	global A, N
	for i in range(N + 1):
		for j in range(N + 1):
			if j < i:
				A[i][j] = -random.random()
			else:
				A[i][j] = 0.0

def set_test_A():
	global A, N
	for i in range(N + 1):
		for j in range(N + 1):
			if j == i - 1:
				A[i][j] = 1.0
			else:
				A[i][j] = 0.0

def set_custom_A():
	global A
	# A = [[0.0, 0.0, 0.0, 0.0], [-0.56345207, 0.0, 0.0, 0.0], [-0.52730513, -0.46009423, 0.0, 0.0], [-0.51084334, -0.3696359, -0.26416908, 0.0]]
	# A = [[0.0, 0.0, 0.0, 0.0], [-0.56347361, 0.0, 0.0, 0.0], [-0.52720981, -0.46007276, 0.0, 0.0], [-0.51079819, -0.36966256, -0.26417685, 0.0]]
	# A = [[0.0, 0.0, 0.0, 0.0, 0.0], [-0.58588916, 0.0, 0.0, 0.0, 0.0], [-0.56640217, -0.53344992, 0.0, 0.0, 0.0], [-0.57019188, -0.46401602, -0.41793027, 0.0, 0.0], [-0.56106712, -0.42894825, -0.31732112, -0.24260608, 0.0]]
	# A = [[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [-0.5956706, 0.0, 0.0, 0.0, 0.0, 0.0], [-0.58295775, -0.56796822, 0.0, 0.0, 0.0, 0.0], [-0.60081277, -0.52712114, -0.50045392, 0.0, 0.0, 0.0], [-0.61200695, -0.50186783, -0.40868126, -0.3843767, 0.0, 0.0], [-0.6029373, -0.48154674, -0.36797587, -0.28231826, -0.23025824, 0.0]]
	A = [[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [-0.59977684, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [-0.58885524, -0.58551022, 0.0, 0.0, 0.0, 0.0, 0.0], [-0.6120996, -0.56242917, -0.54211461, 0.0, 0.0, 0.0, 0.0], [-0.63945231, -0.55535586, -0.48022319, -0.47025945, 0.0, 0.0, 0.0], [-0.64818416, -0.54330915, -0.43877023, -0.36530049, -0.35937924, 0.0, 0.0], [-0.6382574, -0.52817171, -0.41540461, -0.32416287, -0.25897076, -0.22303354,	0.0]]

def minimize():
	global A
	# set_default_A()
	set_custom_A()
	v0 = get_v()

	fun = lambda v: (set_A(v), calc())[-1]

	RES = scipy.optimize.minimize(fun, v0, options={"gtol": 0.00000000001})

	print("=================================================================")
	set_A(RES.x)
	print(calc())
	print(RES)

# set_default_A()
# set_test_A()
# print(A)
# print(calc())
minimize()
# set_test_A()
# calc()
# set_custom_A()
# (a, b) = calc()
# print('%.20f'%a, '%.20f'%b)

# h = 0.000001
# v0 = get_v()
# (y0, y00) = calc()
# der = []
# for i in range(6):
# 	H = 10.0**(-7)
# 	v1 = []
# 	for v in v0:
# 		v1.append(v)
# 	v1[i] += H
# 	set_A(v1)
# 	(y1, y11) = calc()
# 	der.append((y11 - y00)/H)
# 	# print('%.15f'%float((y1 - y0)/H**2), end='\t')
# 	# print("   ", '%.15f'%float((y11 - y00)/H), end='\t')
# 	# print()

# print(der)

# for i in range(15):
# 	H = 10.0**(-i)
# 	v1 = []
# 	for i in range(len(v0)):
# 		v1.append(v0[i] + H * der[i])
# 	# v1[i] += H * der[i]
# 	set_A(v1)
# 	(y1, y11) = calc()
# 	print(y1)

# for i in range(N + 1):
# 	print("i = ", i)
# 	for n in nodes[OFFSET: OFFSET + 2**INT_ORDER + 1]:
# 		print(T*i + n, L2__[i](n) - 2*logcosh(T*i + n))

# for i in range(N + 1):
# 	for j in range(N + 1):
# 		print("i = ", i, "  j = ", j)
# 		for n in nodes[OFFSET: OFFSET + 2**INT_ORDER + 1]:
# 			print(L__0[i][j](n), L__T[i][j](n))



#  N = 1, T = 1.0
#[[ 0.          0.        ]
# [-0.36886773  0.        ]]
# RESULT: 0.48582279393270617


#  N = 2, T = 1.0
# [[ 0.          0.          0.        ]
#  [-0.50939398  0.          0.        ]
#  [-0.44902704 -0.30167512  0.        ]]
# RESULT: 0.30076057861994904


#  N = 3, T = 1.0
# [[ 0.          0.          0.          0.        ]
#  [-0.56345207  0.          0.          0.        ]
#  [-0.52730513 -0.46009423  0.          0.        ]
#  [-0.51084334 -0.3696359  -0.26416908  0.        ]]
# RESULT: 0.19843790699059244
# [[ 0.          0.          0.          0.        ]
#  [-0.56347361  0.          0.          0.        ]
#  [-0.52720981 -0.46007276  0.          0.        ]
#  [-0.51079819 -0.36966256 -0.26417685  0.        ]]



#  N = 4, T = 1.0
# [[ 0.          0.          0.          0.          0.        ]
#  [-0.58588916  0.          0.          0.          0.        ]
#  [-0.56640217 -0.53344992  0.          0.          0.        ]
#  [-0.57019188 -0.46401602 -0.41793027  0.          0.        ]
#  [-0.56106712 -0.42894825 -0.31732112 -0.24260608  0.        ]]
# RESULT: 0.13669968395951987


#  N = 5, T = 1.0
# [[ 0.          0.          0.          0.          0.          0.        ]
#  [-0.5956706   0.          0.          0.          0.          0.        ]
#  [-0.58295775 -0.56796822  0.          0.          0.          0.        ]
#  [-0.60081277 -0.52712114 -0.50045392  0.          0.          0.        ]
#  [-0.61200695 -0.50186783 -0.40868126 -0.3843767   0.          0.        ]
#  [-0.6029373  -0.48154674 -0.36797587 -0.28231826 -0.23025824  0.        ]]
# RESULT: 0.0971609187730787


#  N = 6, T = 1.0
# [[ 0.          0.          0.          0.          0.          0.          0.        ]
#  [-0.59977684  0.          0.          0.          0.          0.			 0.        ]
#  [-0.58885524 -0.58551022  0.          0.          0.          0.			 0.        ]
#  [-0.6120996  -0.56242917 -0.54211461  0.          0.          0.			 0.        ]
#  [-0.63945231 -0.55535586 -0.48022319 -0.47025945  0.          0.			 0.        ]
#  [-0.64818416 -0.54330915 -0.43877023 -0.36530049 -0.35937924  0.			 0.        ]
#  [-0.6382574  -0.52817171 -0.41540461 -0.32416287 -0.25897076 -0.22303354	 0.        ]]
# RESULT: 0.07072350821421666


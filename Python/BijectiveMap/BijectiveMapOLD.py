# kernel is: <f(x0), f(x1)> = exp(-||x0 - x1||^2)
# Let M(x) be the solution of ODE y(0) = x and y'(t) = f(y(t)) and M(x) = y(T)
# Then we wish to calculate M^{-1}(0) using evaluations of M(x) at adaptovely chosen points.

# THIS VERSION IS OLD, UNFINISHED AND INCORRECT

import scipy
import scipy.integrate
import math

N = 1
T = 1.0
# A is strictly lower diagonal
A = [[0.0 for j in range(N + 1)] for i in range(N + 1)]
EPS = 1e-100

# ==================================================================================

# Two key functions
L2 = [lambda t: 0.0 for i in range(N + 1)]
L = [[lambda t: 0.0 for j in range(N + 1)] for i in range(N + 1)]

# Helper functions
L2_at_0 = [0.0 for i in range(N + 1)]
L_at_0 = [[0.0 for j in range(N + 1)] for i in range(N + 1)]
L_at_T = [[0.0 for j in range(N + 1)] for i in range(N + 1)]

V_tilde_at_0 = [[0.0 for j in range(N + 1)] for i in range(N + 1)]


def logcosh(t):
	t = abs(t)
	return t + math.log((1 + math.exp(-2*t))/2.0)

def get_S_integrated(V_tilde, n):
	global A
	S = [0.0 for _ in range(n + 1)]

	for j in range(n + 1):
		S[n][j] = 0.0
		for k in range(j):
			S[n][j] += A[j][k] * (V_tilde[k] - S[n][k])
	
	return S

def get_func(t, y):
	global A
	n = len(y) - 2

	res = [0.0 for _ in range(n + 2)]
	S = [0.0 for _ in range(n + 1)]
 
	for j in range(n + 1):
		S[j] = 0.0
		for k in range(j):
			S[j] += A[j][k] * (S[k] + math.exp(-y[n + 1] + y[k]))
		res[j] = 2*S[j] + math.exp(-y[n + 1] + y[k])
	
	res[n + 1] = 2*(math.tanh(t) + S[n])
	return res

def get_L_func(t, ode, n, j):
	V_tilde = ode.sol(t)
	S = get_S_integrated(V_tilde, n)
	return V_tilde[j] - S[j]

# calculate everything for A.
def calculate():
	global N, T, A

	# Firstly set everything for n = 0.
	L2[0] = lambda t: 2*logcosh(t)
	L[0][0] = lambda t: logcosh(t) + logcosh(T) - logcosh(T - t)

	L2_at_0[0] = 0.0
	L_at_0[0][0] = 0.0
	L_at_T[0][0] = 2*logcosh(T)


	for n in range(N + 1):
		# At n-th step, we assume that L2[i] and L[i][j] are known
		# for all i, j < n.

		# Firstly update the values of L and L2 at 0.
		L2_at_0[n] = 0.0
		for k in range(n):
			for l in range(n):
				L2_at_0[n] += A[n][k] * A[n][l] * L_at_T[k][l]

		for j in range(n):
			L_at_0[n][j] = 0.0
			for k in range(n):
				L_at_0[n][j] += A[n][k] * L_at_T[k][j]

		# Prepare the starting values of V_tilde.
		V_tilde_at_0[n][n] = math.log(math.tanh(T)) + L2_at_0[n]
		for j in range(n):
			func = lambda t: math.exp(-L2[j](t) + 2*sum([A[n][k] * L[j][k](t) for k in range(n)]))
			y = scipy.integrate.quad(func, 0, T)[0]

			V_tilde_at_0[n][j] = math.log(max(y, EPS))

		# Prepare the starting values of y0.
		# values from 0 to n represent V_tilde[n][j](t).
		# value at n+1 represents L2[n](t)
		y0 = [0.0 for _ in range(n + 1)]

		for j in range(n):
			y0[j] = V_tilde_at_0[n][j]
		y0[n + 1] = L2_at_0[n]

		# Propagate the solution and remember its values.
		ode = scipy.integrate.solve_ivp(get_func, (0, T), y0, method='DOP853', dense_output=True)
		
		# Update the function L2[n]
		L2[n] = lambda t: ode.sol(t)[n + 1]

		# Update the functions L2[n] and L[n][j] for 0 <= j <= n
		S0 = get_S_integrated(V_tilde_at_0[n], n)
		ST = get_S_integrated(ode.y, n)

		for j in range(n + 1):
			L_at_T[n][j] = L_at_0[n][j] + (ode.y[j] - ST[j]) - (y0[j] - S0[j])
			L_at_T[j][n] = L_at_T[n][j]
			L[n][j] = lambda t: L_at_0[n][j] + get_L_func(t, ode, n, j) - (y0[j] - S0[j])

		# Finally, update L[i][n] for 0 <= i <= n - 1
		for i in range(n + 1):
			L_at_0[i][n] = 0.0
			for k in range(i):
				L2_at_0[i][n] += A[i][k] * L_at_T[k][n]




	# At last, return the desired value of ||y_N(T)||^2
	return L2[N](T)


# ==================================================================================


def set_default_A():
	global A, N
	for i in range(N + 1):
		for j in range(N + 1):
			if j < i:
				A[i][j] = -1.0
			else:
				A[i][j] = 0.0


RES = scipy.integrate.solve_ivp((lambda t, y: y), (0.0, 1.0), [1.0], method='DOP853', dense_output=True)
print(RES)

for i in range(11):
	print(math.exp(i/10.0) - RES.sol(i/10.0))

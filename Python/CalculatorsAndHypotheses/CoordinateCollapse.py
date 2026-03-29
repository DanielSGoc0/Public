# Check if a particular kind of coordinate
# collapse is always helpful.
import numpy as np
import scipy
import random
import scipy.linalg

np.set_printoptions(suppress=True, formatter={'float_kind':'{:5.8f} \t'.format}, linewidth=200000, threshold=np.inf)

n = None
m = None
N = None
epsilons = None


# Here is how indices work:
# 0:n corresponds to initial data
# n:m corresponds to first parallel evaluation
# m:N corresponds to latter evaluations
# N is the output point

def phi(r2):
	# return np.exp(-r2 / 2.0)
	return np.exp(-r2 / 4.0) + np.exp(-r2)

# multiplied by -2.0 in fact
def dphi(r2):
	# return np.exp(-r2 / 2.0)
	return np.exp(-r2 / 4.0) / 2.0 + np.exp(-r2) * 2.0

def calc(X, sigma2):
	global n, m, N, epsilons

	Sigma = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N + 1):
		for j in range(N + 1):
			Sigma[i][j] = phi(np.dot(X[i] - X[j], X[i] - X[j]))

	for i in range(N):
		Sigma[i][i] += sigma2[i]

	ANS1 = Sigma[N, :N] @ np.linalg.solve(Sigma[:N, :N], Sigma[:N, N])


	ANS2 = 0.0
	for epsilon in epsilons:
		Sigma = np.zeros((N + 1, N + 1), dtype='double')
		for i in range(N + 1):
			for j in range(N + 1):
				r2_1 = np.dot(X[i][:n] - X[j][:n], X[i][:n] - X[j][:n])
				r2_3 = np.dot(X[i][m:N] - X[j][m:N], X[i][m:N] - X[j][m:N])
				r_i = np.sqrt(np.dot(X[i][n:m], X[i][n:m]))
				r_j = np.sqrt(np.dot(X[j][n:m], X[j][n:m]))

				Sigma[i][j] = phi(r2_1 + (epsilon[i] * r_i - epsilon[j] * r_j)**2 + r2_3)

		for i in range(N):
			Sigma[i][i] += sigma2[i]

		ANS2 = max(ANS2, Sigma[N, :N] @ np.linalg.solve(Sigma[:N, :N], Sigma[:N, N]))
		

	return (ANS1, ANS2)


def full_derivative(X, sigma2):
	global m, N

	Sigma = np.zeros((N + 1, N + 1), dtype='double')
	dSigma = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N + 1):
		for j in range(N + 1):
			Sigma[i][j] = phi(np.dot(X[i] - X[j], X[i] - X[j]))
			dSigma[i][j] = dphi(np.dot(X[i] - X[j], X[i] - X[j]))

	for i in range(N):
		Sigma[i][i] += sigma2[i]

	b = np.linalg.solve(Sigma[:N, :N], Sigma[:N, N])

	A = np.copy(dSigma[m:N, m:N])
	for k in range(m, N):
		A[k - m, k - m] += (dSigma[k][N] - dSigma[k, :N] @ b) / b[k]

	print(A)
	return np.linalg.eigvalsh(A) 



n = 2
m = 4
N = 6

epsilons = []
e = [1 for i in range(N + 1)]
while True:
	epsilons.append(e)

	END = True
	for i in range(m, N):
		if e[i] == 1:
			e[i] = -1
			END = False
			break
		e[i] = 1
	if END:
		break


for tries in range(1000000):
# for tries in range(1):

	X = np.zeros((N + 1, N + 1), dtype='double')
	sigma2 = np.zeros((N), dtype='double')

	for i in range(N + 1):
		for j in range(n):
			X[i][j] = random.normalvariate(0, 1.0)
	for i in range(m, N + 1):
		for j in range(n, m):
			X[i][j] = random.normalvariate(0, 1.0)
	for j in range(m, N):
		X[N][j] = random.normalvariate(0, 1.0)
	
	for i in range(N):
		# sigma2[i] = 0.0
		sigma2[i] = random.random() / (1.0 - random.random())
		# sigma2[i] = 1.0
	# sigma2[0] = 1.0

	# (ANS1, ANS2) = calc(X, sigma2)

	# print(tries)
	# if ANS1 > ANS2:
	# 	print(ANS1, ANS2)
	# 	print("dang")
	# 	print(X)
	# 	print(sigma2)
	# 	break
	
	evals = full_derivative(X, sigma2)
	
	WORKS = True
	for eval in evals:
		if eval < 0.0:
			WORKS = False
			break
	
	if not WORKS:
		print("dang")
		print(evals)
		print(X)
		print(sigma2)
		break

	print(tries, '\t', evals)
	# print(X)



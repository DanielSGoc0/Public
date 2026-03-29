# Calculates the "effectiveness" of infinite-size "partially" fixed
# coordinates method, especially in the errored setting

import numpy as np
import scipy
import random
import pprint

import scipy.linalg

np.set_printoptions(suppress=True,
   formatter={'float_kind':'{:5.8f}\t'.format}, linewidth=20000)

err2 = 0.0
COEFF = 1.0
DEBUG = False
M = 1024

def get_v(N, x, x0, c, c0):

	v = []
	for i in range(1, N + 1):
		v.append(x[i])
	v.append(x0)
	for i in range(N):
		v.append(c[i])
	v.append(c0)

	return v


def get_v2(N, x):

	v = []
	for i in range(1, N + 1):
		v.append(x[i])

	return v


def get_x(N, v):
	x = [0.0]
	c = []
	for i in range(N):
		x.append(v[i])
	x0 = v[N]
	for i in range(N):
		c.append(v[i + N + 1])
	c0 = v[2*N + 1]

	return (x, x0, c, c0)


def get_x2(N, v):
	x = [0.0]
	for i in range(N):
		x.append(v[i])

	return x


# ===============================================================================

def sqr(x):
	if x > 1:
		return x
	elif x < -1:
		return -x
	else:
		return x**2

# calculates log(1 + (e^(-d) - i)^2) + log(1 + (1 - e^d * i)^2)
def error(d, i):
	if i > 0:
		return (np.exp(-d/2.0) - i)**2 + (-d/2.0 - np.log(i))**2
	else:
		return i**2

# X is strictly lower diagonal
def calc(t):
	global err2, COEFF, DEBUG

	if DEBUG:
		print("--------------------")

	x = np.array(t[0])
	x0 = t[1]
	c = np.array(t[2])
	c0 = t[3]
	N = len(t[0]) - 1
	q = 1.0/(1.0 - np.exp(-x0*x0))

	if DEBUG:
		print(x)
		print(x0)
		print(c)
		print(c0)

	dists = np.zeros((N + 1), dtype='double')
	inners = np.zeros((N + 1), dtype='double')

	for s in range(N + 1):
		for i in range(-N, -N + s):
			dists[s] += (x[-i] - x0)**2
		for i in range(-N + s, 0):
			dists[s] += (x[-i] - x[-i + s])**2
		for i in range(0, s):
			dists[s] += x[-i + s]**2

	for s in range(N + 1):
		inners[s] = np.exp(-x0*x0 * s / 2.0) * c0*c0 * q
		for k in range(-N + 1, -N + s + 1):
			inners[s] += c0 * c[-k] * np.exp(-(s - N - k) * x0*x0 / 2.0)
		for k in range(-N + s + 1, 1):
			inners[s] += c[-k] * c[-k + s]
	inners[0] -= err2
	
	# RES = x0 * c0 / (np.exp(x0*x0 / 2.0) - 1.0)
	RES = 0.0
	sum_c = 0.0
	for k in range(N):
		sum_c += c[k]**2 
	if sum_c < 1.0 + err2:
		RES = np.exp(-x0*x0 / 2.0) * np.sqrt((np.exp(x0*x0/2.0) + 1) * (x0 * x0 / (np.exp(x0*x0 / 2.0) - 1.0)) * (1.0 + err2 - sum_c))

	RES += x[N] * c0
	for k in range(-N + 1, 1):
		RES += x[-k] * c[-k]
	RES *= 2.0*(1.0 + err2) / (sum_c + c0*c0 * q + 1.0 + err2)
	if DEBUG:
		print("modulus:", sum_c + c0*c0 * q)
		print("RES1:", RES)
	
	if sum_c > 1.0 + err2:
		RES -= COEFF * (1.0 + err2 - sum_c)**2

	# Our goal is to maximize RES under the constraints inners[s] = exp(-dists[s] / 2) + kronecker(0, s) * error2 for all 0 <= s <= N

	sum_x = 0.0
	for k in range(N + 1):
		sum_x += x[k]**2
	sum_x += x0**2
	for s in range(N + 1):
		diff2 = (np.exp(-dists[s]/2) - inners[s])**2
		if DEBUG:
			print("s =", s, dists[s], inners[s], diff2)
		# RES -= COEFF * (np.log(1 + diff2) + np.log(np.exp(-dists[s]) + diff2) + dists[s])
		RES -= COEFF * diff2 * (1.0 + np.exp(dists[s])) * (1 + sum_x)
		# RES -= COEFF * error(dists[s], inners[s])
	if DEBUG:
		print("RES2:", RES)
	return RES

# X is strictly lower diagonal
def calc2(t):
	global err2, COEFF, DEBUG, M

	if DEBUG:
		print("--------------------")

	x = np.array(t)
	print(x)
	N = len(t) - 1
	
	X = np.zeros((M, M), dtype='double')

	if False:
		for i in range(M):
			X[i][i] = 1.0 + err2

			d = 0.0
			for k in range(i):
				d += x[min(k, N)]**2
			X[i][0] = np.exp(-d / 2.0)
			X[0][i] = X[i][0]

			for j in range(1, i):
				X[i][j] = X[i - 1][j - 1] * np.exp(-(x[min(i, N)] - x[min(j, N)])**2 / 2.0)
				X[j][i] = X[i][j]

	else:

		dists = np.zeros((N), dtype='double')

		for s in range(N):
			for i in range(-N + 1, -N + 1 + s):
				dists[s] += (x[-i] - x[N])**2
			for i in range(-N + 1 + s, 0):
				dists[s] += (x[-i] - x[-i + s])**2
			for i in range(0, s):
				dists[s] += x[-i + s]**2

		for i in range(M):
			for j in range(M):
				d = abs(i - j)
				if d == 0:
					X[i][j] = 1.0 + err2
				elif d < N:
					X[i][j] = np.exp(-dists[d] / 2.0)
				else:
					X[i][j] = np.exp(-dists[N - 1] / 2.0 - (d - N + 1) * x[N]**2 / 2.0)

	
	C = scipy.linalg.cholesky(X, lower=True)
	print(C)
	print(np.linalg.inv(C)[M - 1])
	
	RES = 0.0
	for k in range(min(N, M)):
		RES += C[M - 1][M - 1 - k] * x[k]
	for k in range(min(N, M), M):
		RES += C[M - 1][M - 1 - k] * x[N]
	return RES

# M is the DFT size
def calc3(t):
	global err2, COEFF, DEBUG, M

	x = np.array(t)
	N = len(t) - 1
	# x[1] = 1.0
	# x[N] = max(x[N], 0.1)
	x[N] = 1.0

	dists = np.zeros((N), dtype='double')

	for s in range(1, N):
		for i in range(-N + 1, -N + 1 + s):
			dists[s] += (x[-i] - x[N])**2
		for i in range(-N + 1 + s, 0):
			dists[s] += (x[-i] - x[-i + s])**2
		for i in range(0, s):
			dists[s] += x[-i + s]**2

	# We also find dists[N - 1 + s] = dists[N - 1] + s * x[N]**2

	# Now, we begin by finding the symbol:

	symbol = np.zeros(M, dtype='complex')
	symbol[0] = 1.0 + err2
	for s in range(1, N - 1):
		v = np.exp(-dists[s] / 2.0)
		symbol[s % M] += v
		symbol[-(s % M)] += v
	for a in range(M):
		v = np.exp(-(dists[N - 1] + a * x[N]**2) / 2.0) / (1.0 - np.exp(-M * x[N]**2 / 2.0))
		symbol[(a + N - 1) % M] += v
		symbol[-((a + N - 1) % M)] += v
	symbol = np.real(np.fft.fft(symbol))

	# With symbol calculated, we apply the Hilbert transformation

	hilbert_transform = np.zeros(M, dtype='complex')
	for k in range(M):
		if symbol[k] < 0.00000000001:
			print(">:(")
			hilbert_transform[k] = -12.664218011
		else:
			hilbert_transform[k] = np.log(symbol[k]) / 2.0
	hilbert_transform = np.fft.fft(hilbert_transform)

	hilbert_transform[0] = 0.0
	hilbert_transform[M//2] = 0.0
	for k in range(1, M//2):
		hilbert_transform[k] *= -1.0j
		hilbert_transform[-k] *= 1.0j

	# The Fourier Transform of Cholesky decomposition is a direct function of hilbert_transform and the symbol

	c_hat = np.fft.ifft(hilbert_transform)
	for k in range(M):
		c_hat[k] = np.exp(np.real(c_hat[k]) * 1.0j)
		if symbol[k] < 0.00000000001:
			c_hat[k] = 0.0
		else:
			c_hat[k] *= np.sqrt(symbol[k])
	
	# Finally, we wish to compute the sum c[k] * x[k]

	c = np.real(np.fft.ifft(c_hat))
	c[1:] = c[1:][::-1]
	print(c[:20])
	RES = 0.0
	for k in range(1, N):
		RES += x[k] * c[k % M]
	for k in range(N, M):
		RES += x[N] * c[k]

	print(x, RES)

	return RES



# =========================================================================================

def initial_guess(N):
	x = [0.0]
	x0 = 1.0
	c = []
	c0 = 1.0
	for i in range(N):
		x.append(random.random())
		c.append(random.random())
	return (x, x0, c, c0)

def initial_guess2(N):
	x = [0.0]
	for i in range(N):
		# x.append(random.random())
		x.append(i/N)
	return x


def minimize():
	global err2, COEFF, DEBUG, M
	N = 100
	M = 32768
	# (x, x0, c, c0) = initial_guess(N)
	x = initial_guess2(N)
	# x = [0.0, 0.18687467386503206, 0.18765184710700153, 0.1885431350999313, 0.18951559058692458, 0.19054323324570027, 0.19160555681627264, 0.19268727221783255, 0.19377649960566223, 0.19486453795250583, 0.1959450448357043, 0.19701369453787865, 0.1980679060616212, 0.19910638182850782, 0.20012932722771928, 0.20113817745486406, 0.20213551891070752, 0.2031251731697674, 0.2041121857743966, 0.20510306472439216, 0.21560086021172217]
	# x0 = 1.0
	# x = [0.0, 1.0, 1.0]
	# c0 = 1.0
	# c = [1.0, 1.0]


	fun = lambda v: -calc3(get_x2(N, v))
	COEFF = 1.0

	# for i in range(10):
	# 	print("=========================   ", i, "   ==========================")
	# 	RES = scipy.optimize.minimize(fun, get_v(N, x, x0, c, c0), method='Nelder-Mead')
	# 	t = get_x(N, RES.x)
	# 	x = t[0]
	# 	x0 = t[1]
	# 	c = t[2]
	# 	c0 = t[3]
	# 	print(x, x0, c, c0)
	# 	print(RES)

	# 	DEBUG = True
	# 	fun(get_v(N, x, x0, c, c0))
	# 	DEBUG = False

	# 	COEFF *= 2


	err2 = 1.0
	RES = scipy.optimize.minimize(fun, get_v2(N, x), options={"gtol": 0.0000000000000000001})
	print(RES)
	x = get_x2(N, RES.x)
	print(x)

minimize()

# N = 2
# (x, x0, c, c0) = initial_guess(N)
# calc2((x, x0, c, c0))

# M = 64

# RES = 0.0
# for k in range(M):
# 	x = (2.0*k + 1.0) / (2.0 * M)
# 	RES += np.log(1.0 - np.cos(x)) * np.cos((1.0 - x) / 2.0) / np.sin((1.0 - x) / 2.0)
# print(RES / M)

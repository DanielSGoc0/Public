# Used to verify an old hypothesis: can convolving by a non-negative
# distribution by a Gaussian kernel leave it somehow "invariant"?

# turns out the hypothesis is wrong!

import math
import scipy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator

N = 24
M = 100
S = 180
gamma = -0.25
LIMIT = min(S, math.floor(-1.0/gamma) + 1)
rgammas = [scipy.special.rgamma(-gamma*k) for k in range(S)]

h = 0.01

factor = [1.0 for i in range(S)]
term = [1.0 for i in range(S)]
# factor[k] = c[k] / c[k - 1], with the assumption c[-1] = 1
# factor[0] * factor[1] * ... factor[k] * term[k] = c[k] * term[k]
# the rule is: factor[0] * (term[0] + a * factor[1] * (term[1] + a * factor[2] * (term[2] + a * 0)))


tab = [[0.0 for j in range(-M, M + 1)] for i in range(-N, N + 1)]


# Gaussian distribution
# for k in range(S//2):
# 	if k == 0:
# 		factor[2*k] = 1.0
# 		factor[2*k + 1] = 1.0
# 	else:
# 		factor[2*k] = 1.0/k
# 		factor[2*k + 1] = 1.0
# 	term[2*k] = 1.0
# 	term[2*k + 1] = 0.0
# for k in range(S - 2):
# 	term[k] = 0.0

# Laplace distribution
# factor[0] = 1
# factor[1] = 1
# factor[2] = 1
# factor[3] = 1
# factor[4] = 9
# factor[5] = 1
# term[0] = 1
# term[1] = 0
# term[2] = 2
# term[3] = 0
# term[4] = 1
# term[5] = 0
# for k in range(3, S//2):
# 	factor[2*k] = 4.0
# 	factor[2*k + 1] = 1.0
# 	term[2*k] = 1.0
# 	term[2*k + 1] = 0.0

# log logistic
factor[0] = 8/math.pi**2
for k in range(1, S):
	factor[k] = 2
for k in range(0, S//2):
	term[2*k] = (1 - 2**(-k - 2)) * scipy.special.zeta(k + 2)
	term[2*k + 1] = 0.0
# for k in range(S - 2):
# 	term[k] = 0.0

# beta-prime with 1/2, 1/2
factor[0] = 1.0/2
term[0] = 2
term[1] = 0
for k in range(1, S):
	factor[k] = (2*math.pi)/k
lis = scipy.special.euler(S)
for k in range(1, S//2):
	term[2*k] = (-1)**k * lis[2*k]
	term[2*k + 1] = 0.0
# for k in range(S - 2):
# 	term[k] = 0.0


def integrand(a, z, s, L, S):
	global rgammas, factor, term, gamma

	w = a * s**(-gamma)
	inner = 0.0
	for k in range(S - 1, L - 1, -1):
		inner *= w
		inner += term[k] * rgammas[k]
		inner *= factor[k]
	inner *= s**(-gamma*L - 1)
	inner *= math.exp(-s)
	inner *= math.cos(z * math.sqrt(2*s))

	return inner

def calculate(a, z, L, S):
	global N, M, gamma, factor, term

	RES = 0.0
	for k in range(L - 1, -1, -1):
		# print(TOTAL)
		RES *= a
		RES += term[k] * scipy.special.hyp1f1(-gamma*k, 1.0/2, -z*z/2.0)
		RES *= factor[k]
	
	integral = scipy.integrate.quad(lambda s: integrand(a, z, s, L, S), 0, np.inf)
	# print(integral)
	integral = integral[0]
	for k in range(L - 1, -1, -1):
		integral *= factor[k] * a

	RES += integral

	return RES


print(calculate(0.238, 0.8, LIMIT, S))

# for i in range(-N, N + 1):
# 	for j in range(-M, M + 1):
# 		a = i*h
# 		z = j*h

# 		tab[i + N][j + M] = calculate(a, z, LIMIT, S)

# for i in range(1, 2*N):
# 	for j in range(1, 2*M):
# 		a = (i - N)*h
# 		z = (j - M)*h

# 		print((tab[i][j + 1] + tab[i][j - 1] - 2*tab[i][j])/(h*h) + z * (tab[i][j + 1] - tab[i][j - 1])/(2*h) - 2*a*gamma*(tab[i + 1][j] - tab[i - 1][j])/(2*h))


# for i in range(0, 2*N + 1):
# 	a = (i - N)*h

# 	# print(tab[i][M] - (1 - a*a)**2 / (1 - 4*a**2))
# 	# if a == 0:
# 	# 	# print(tab[i][M] - 1)
# 	# 	print(tab[i][M], 1)
# 	# else:
# 	# 	print(tab[i][M], math.tan(a*math.pi)/(a*math.pi))
# 	# 	# print(tab[i][M] - math.tan(a*math.pi)/(a*math.pi)
# 	print(tab[i][M], 1.0/2*(1 + 1.0/math.cos(2*math.pi*a)))
# 	# print(tab[i][M] - math.exp(a*a))

# fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
# ax.autoscale(enable=None, axis="x", tight=True)

# X = [[h * j for j in range(-M, M + 1)] for i in range(-N, N + 1)]
# Y = [[h * i for j in range(-M, M + 1)] for i in range(-N, N + 1)]
# Z = np.array(tab)

# surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm, linewidth=0, antialiased=False)

# plt.show()

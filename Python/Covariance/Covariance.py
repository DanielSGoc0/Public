# Used to verify an old hypothesis: can convolving by a non-negative
# distribution by a Gaussian kernel leave it somehow "invariant"?

# turns out the hypothesis is wrong!

import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator

N = 100
M = 100
gamma = -0.25

s = 0
t = 3
h = (t - s)/(N + 2*M)

def A(a):
	return math.exp(a*a/2)
	# return 1/(1 - a*a)
	# return 3*a*a + 14*a + 7

def quot(a):
	return a
	# return 2*a/(1 - a*a)
	# return (6*a + 14)/(3*a*a + 14*a + 7)

tab = [[0.0 for j in range(M + 1)] for i in range(N + 2*M + 1)]

for i in range(N + 2*M + 1):
	a = s + i*h
	tab[i][0] = A(2*a)

for j in range(M):
	for i in range(j + 1, N + 2*M + 1 - j - 1):
		a = s + h*i
		z = h*j
		
		# (tab[i][j + 1] + tab[i][j - 1] - 2*tab[i][j])/h**2/2 = -2*a*gamma*quot(a)*tab[i][j] + a*gamma*(tab[i + 1][j] - tab[i - 1][j])/(2*h) - z/2.0 * (tab[i][j + 1] - tab[i][j - 1])/(2*h)
		# tab[i][j + 1] * (1 + z*h/2) = tab[i][j]*2*(1 - 2*a*gamma*quot(a)*h**2) + a*gamma*(tab[i + 1][j] - tab[i - 1][j])*h - tab[i][j - 1] * (1 - z*h/2)

		# tab[i][j + 1] * 2 = tab[i][j]*2*(1 - 2*a*gamma*quot(a)*h**2) + a*gamma*(tab[i + 1][j] - tab[i - 1][j])*h
		# tab[i][j + 1] * 2 = tab[i][j]*2*(1 - a*gamma*quot(a)*h**2) + a*gamma*(tab[i + 1][j] - tab[i - 1][j])*h/2.0


		if j == 0:
			tab[i][j + 1] = tab[i][j]*2*(1 - 2*a*gamma*quot(a)*h**2) + a*gamma*(tab[i + 1][j] - tab[i - 1][j])*h
			tab[i][j + 1] /= 2
		else:
			tab[i][j + 1] = tab[i][j]*2*(1 - 2*a*gamma*quot(a)*h**2) + a*gamma*(tab[i + 1][j] - tab[i - 1][j])*h - tab[i][j - 1] * (1 - z*h/2)
			tab[i][j + 1] /= (1 + z*h/2.0)
			# print((tab[i][j + 1] + tab[i][j - 1] - 2*tab[i][j])/h**2 - (-2*a*gamma*quot(a)*tab[i][j] + a*gamma*(tab[i + 1][j] - tab[i - 1][j])/(2*h) - z/2.0 * (tab[i][j + 1] - tab[i][j - 1])/(2*h)))






fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
ax.autoscale(enable=None, axis="x", tight=True)

X = [[h * j for j in range(M + 1)] for i in range(M, N + M + 1)]
Y = [[h * i for j in range(M + 1)] for i in range(M, N + M + 1)]
Z = np.array([[tab[i][j] for j in range(M + 1)] for i in range(M, N + M + 1)])

surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm, linewidth=0, antialiased=False)

plt.show()

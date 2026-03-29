# Original 2D plotter

import itertools
import scipy
import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
from matplotlib.widgets import Button


# ==============================================================================================================================================


def area(v0, v1, v2):
	return np.dot(v0 - v1, v0 - v1) * np.dot(v0 - v2, v0 - v2) - np.dot(v0 - v1, v0 - v2)**2

def custom_sum(X, s):
	N = X.shape[0]

	S = 0.0
	for perm in itertools.permutations(s):
		Sigma = np.zeros((N, N), dtype = 'double')

		for i in range(N):
			for j in range(N):
				Sigma[i][j] = (-np.dot(X[i] - X[j], X[i] - X[j])/2)**perm[i] / scipy.special.factorial(perm[i])

		S += np.linalg.det(Sigma)
		
	return S

def F(x, y):
	# N = 4
	# D = 2

	# X = np.zeros((N, D), dtype='double')
	# X[0] = [0, 0]
	# X[1] = [1, 0]
	# X[2] = [0, 1]
	# X[3] = [x, y]

	# w0 = custom_sum(X, [0, 1, 1, 2]) / 2
	# w1 = custom_sum(X, [1, 1, 1, 1]) / 24
	
	# return np.sqrt(w0 + w1)

	# X = np.array([[x, y], [1.06, 0], [0.06, 0]])
	# Sigma = np.zeros((3, 3), dtype='double')
	# for i in range(3):
	# 	for j in range(3):
	# 		Sigma[i][j] = np.exp(-np.dot(X[i] - X[j], X[i] - X[j]) / 2.0)

	# return ((Sigma[0][2] - Sigma[0][1] * Sigma[1][2])/np.sqrt(1.0 - Sigma[0][1]**2))**2

	return (np.exp(-(x-y)**2 / 2.0) - np.exp(-x**2 / 2.0 - y**2 / 2.0) * (1 - 0.1)) / (1 - np.exp(-x**2) * (1 - 0.1))

# ========================================================================================================================

MAX = 400
MIN_X = -3
MAX_X = 3
MIN_Y = -3
MAX_Y = 3

X = [[MIN_X + (MAX_X - MIN_X) * i / MAX for j in range(MAX + 1)] for i in range(MAX + 1)]
Y = [[MIN_Y + (MAX_Y - MIN_Y) * j / MAX for j in range(MAX + 1)] for i in range(MAX + 1)]
Z = [[0.0 for j in range(MAX + 1)] for i in range(MAX + 1)]

def update(event):
	global X, Y, Z, MAX
	# whatever the button does
	print("magic")


	Z = np.array([[F(X[i][j], Y[i][j]) for j in range(MAX + 1)] for i in range(MAX + 1)])
	ax.clear()
	ax.plot_surface(X, Y, Z, cmap=cm.coolwarm, linewidth=0, antialiased=False)
	plt.draw()

# ========================================================================================================================



fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
ax.autoscale(enable=None, axis="x", tight=True)

Z = np.array([[F(X[i][j], Y[i][j]) for j in range(MAX + 1)] for i in range(MAX + 1)])

surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm, linewidth=0, antialiased=False)

axbutton = fig.add_axes([0.81, 0.05, 0.1, 0.075])
button = Button(axbutton, 'magic')
button.on_clicked(update)

ax.set_xlabel('x')
ax.set_ylabel('y')

plt.show()

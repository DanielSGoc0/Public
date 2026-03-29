# Plotter2D with sliders

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.widgets import Button, Slider

A = 0.25
R = 1.0

# finds an approximate solution of f(x) = y
# f is assumed increasing and we use N iterations
# We search in interval (a0, b0) and so we assume f(a0) <= v <= f(b0)
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

def f(x, y):
	global A, R

	x *= R
	y *= np.sqrt(A)
	
	A0 = max(0.0, 1.0 - np.exp(-x**2) * (1.0 - A + y**2))
	R0 = np.sqrt(R**2 - x**2)
	# A0 = A
	# R0 = R
	# A0 = y
	# R0 = x

	func = lambda z: np.sqrt(4.0*z*np.exp(z) - 4.0*z)/(np.sqrt(A0) + np.sqrt(A0 + 4.0*z*np.exp(z))) - np.sqrt(R0**2 - z)
	z0 = binsearch(func, 0.0, 0.0, R0**2, 50)
	v = np.exp(-R0**2 / 2.0) * (np.sqrt(z0 * A0) + np.sqrt((R0**2 - z0) * (np.exp(z0) - 1)))

	return (x, y, np.exp(-R**2 / 2.0)*x*y + v)
	# return (x, y, v)
	# return (x, y, v - np.exp(-R**2 / 2.0) * (F0 + x*y + np.sqrt((R**2 - x**2) * (np.exp(x**2) - (1.0 - A + y**2)))))
	# return (v - np.exp(-R**2/2.0) * (np.exp(-x**2 / 2.0) * (F + x*y)), 1.0 - np.exp(-x**2)*(1.0 - A + y**2), R**2 - x**2)

# ========================================================================================================================

MAX = 10
MIN_X = 0.0
MAX_X = 1.0
MIN_Y = 0.0
MAX_Y = 1.0

X = [[MIN_X + (MAX_X - MIN_X) * i / MAX for j in range(MAX + 1)] for i in range(MAX + 1)]
Y = [[MIN_Y + (MAX_Y - MIN_Y) * j / MAX for j in range(MAX + 1)] for i in range(MAX + 1)]
# Z = [[0.0 for j in range(MAX + 1)] for i in range(MAX + 1)]

def update(event):
	global X, Y, Z, MAX, ax, A, R, A_slider, R_slider, MIN_X, MAX_X, MIN_Y, MAX_Y

	A = A_slider.val
	R = R_slider.val

	print("EVENT")

	Z = np.array([[(X[i][j], Y[i][j], 0.0) for j in range(MAX + 1)] for i in range(MAX + 1)])

	# for a in range(MAX + 1):
	# 	for b in range(MAX + 1):
	# 		# R = (MIN_X + (MAX_X - MIN_X) * (2.0*a + 1.0) / (2.0*MAX))*0.1
	# 		# A = MIN_Y + (MAX_Y - MIN_Y) * (2.0*b + 1.0) / (2.0*MAX)

	# 		W = np.array([[f(np.sqrt(X[i][j]), Y[i][j]) for j in range(MAX + 1)] for i in range(MAX + 1)])

	# 		supremum = W[0][0][2]
	# 		supremum_at = (0, 0)
	# 		for i in range(MAX + 1):
	# 			for j in range(MAX + 1):
	# 				if W[i][j][2] > supremum:
	# 					supremum = W[i][j][2]
	# 					supremum_at = (i, j)

	# 		Z[a][b] = (R, A, Y[supremum_at[0]][supremum_at[1]])


	Z = np.array([[f(X[i][j], Y[i][j]) for j in range(MAX + 1)] for i in range(MAX + 1)])
	ax.clear()
	ax.set_xlabel('R Label')
	ax.set_ylabel('A Label')
	ax.set_zlabel('value Label')
	ax.plot_surface(Z[:, :, 0], Z[:, :, 1], Z[:, :, 2], cmap=cm.coolwarm, linewidth=0, antialiased=False)
	plt.draw()

	supremum = Z[0][0][2]
	supremum_at = (Z[0][0][0], Z[0][0][1])
	for i in range(MAX + 1):
		for j in range(MAX + 1):
			if Z[i][j][2] > supremum:
				supremum = Z[i][j][2]
				supremum_at = (X[i][j], Y[i][j])

	print("supremum:", supremum)
	print("supremum at:", supremum_at)

# ========================================================================================================================



fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
ax.autoscale(enable=None, axis="x", tight=True)
ax.set_xlabel('R Label')
ax.set_ylabel('A Label')
ax.set_zlabel('value Label')

Z = np.array([[f(X[i][j], Y[i][j]) for j in range(MAX + 1)] for i in range(MAX + 1)])

surf = ax.plot_surface(Z[:, :, 0], Z[:, :, 1], Z[:, :, 2], cmap=cm.coolwarm, linewidth=0, antialiased=False)

axbutton = fig.add_axes([0.81, 0.05, 0.1, 0.075])
button = Button(axbutton, 'magic')
button.on_clicked(update)

# Make a vertical slider to control variable A
axis_A = fig.add_axes([0.1, 0.25, 0.0225, 0.63])
A_slider = Slider(
    ax=axis_A,
    label='A',
    valmin=0.0,
    valmax=1.0,
    valinit=A,
    orientation="vertical"
)

# Make a vertical slider to control variable R
axis_R = fig.add_axes([0.15, 0.25, 0.0225, 0.63])
R_slider = Slider(
    ax=axis_R,
    label='R',
    valmin=0.0,
    valmax=3.0,
    valinit=R,
    orientation="vertical"
)

A_slider.on_changed(update)
R_slider.on_changed(update)

plt.show()

# Plotter2D with sliders

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.widgets import Button, Slider

F = 0.0
A = 1.0

def f2(x, y):
	global F, A

	# y *= np.sqrt(A)

	# F_new = np.exp(-x*x / 2.0) * (F + x * y)
	# w2_new = 1.0 - np.exp(-x*x) * (1.0 - A + y**2)
	F_new = x
	w2_new = y

	x_opt = (np.sqrt(F_new**2 + 4.0 * w2_new) - F_new)/(2.0 * np.sqrt(w2_new))
	y_opt = np.sqrt(w2_new)
	return np.exp(-x_opt*x_opt / 2.0) * (F_new + x_opt*y_opt)

def f(R2, A):
	# X = np.array([[0.0, 0.0], [x, y], [1.0, 0.0]])
	# M = np.zeros((3, 3), dtype='double')
	# for i in range(3):
	# 	for j in range(3):
	# 		M[i][j] = np.exp(-np.dot(X[i] - X[j], X[i] - X[j])/2.0)
	# # M[1][1] += 0.1

	# return np.linalg.det(M) / np.linalg.det(M[0:2, 0:2])
	# # return (f2(x, y + 0.001) - f2(x, y))/(f2(x + 0.001, y) - f2(x, y))
	return np.sqrt(R2 * (1.0 - np.exp(-R2) * (1.0 - A)))

# ========================================================================================================================

MAX = 100
MIN_X = 0.0
MAX_X = 100.0
MIN_Y = 0.0
MAX_Y = 1.0

X = [[MIN_X + (MAX_X - MIN_X) * i / MAX for j in range(MAX + 1)] for i in range(MAX + 1)]
Y = [[MIN_Y + (MAX_Y - MIN_Y) * j / MAX for j in range(MAX + 1)] for i in range(MAX + 1)]
Z = [[0.0 for j in range(MAX + 1)] for i in range(MAX + 1)]

def update(event):
	global X, Y, Z, MAX, ax, F, A, F_slider, A_slider

	F = F_slider.val
	A = A_slider.val

	# whatever the button does
	print("magic")


	Z = np.array([[f(X[i][j], Y[i][j]) for j in range(MAX + 1)] for i in range(MAX + 1)])
	ax.clear()
	ax.plot_surface(X, Y, Z, cmap=cm.coolwarm, linewidth=0, antialiased=False)
	plt.draw()

# ========================================================================================================================



fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
ax.autoscale(enable=None, axis="x", tight=True)

Z = np.array([[f(X[i][j], Y[i][j]) for j in range(MAX + 1)] for i in range(MAX + 1)])

surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm, linewidth=0, antialiased=False)

axbutton = fig.add_axes([0.81, 0.05, 0.1, 0.075])
button = Button(axbutton, 'magic')
button.on_clicked(update)

# Make a vertical slider to control variable F
axis_F = fig.add_axes([0.05, 0.25, 0.0225, 0.63])
F_slider = Slider(
    ax=axis_F,
    label='F',
    valmin=-3.0,
    valmax=3.0,
    valinit=F,
    orientation="vertical"
)

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

F_slider.on_changed(update)
A_slider.on_changed(update)

plt.show()

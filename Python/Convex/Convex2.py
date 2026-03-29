# A program that tries to iteratively solve the Monge-Ampere equation.
# With 2D plotting

import math
import scipy
import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
from matplotlib.widgets import Button

S3 = math.sqrt(3)/2

def F(x, y):
	# return (2*x**2 + 3*y**2 - 2*x*y)/2.0 + 0.529729 + 0.9472834982 * x - 0.354235239 * y
	return math.sqrt((x - 4)**2 + y**2) + math.sqrt((x + 4)**2 + y**2)

def det(x, y, h):
	global S3

	w = F(x, y)
	a_d = F(x, y + h) + F(x, y - h)
	b_e = F(x + S3*h, y + h/2) + F(x - S3*h, y - h/2)
	c_f = F(x + S3*h, y - h/2) + F(x - S3*h, y + h/2)

	return ((a_d + b_e + c_f - 2*w)**2 + 8*w**2 - 2*(a_d**2 + b_e**2 + c_f**2))/(3.0 * h**4)

# ==============================================================================================================================================


N = 20
MAX = 2*N + 1
h = 0.5

X = [[S3 * h * i for j in range(MAX)] for i in range(MAX)]
Y = [[(j - i/2.0) * h for j in range(MAX)] for i in range(MAX)]
Z = [[0.0 for j in range(MAX)] for i in range(MAX)]
# X = [[S3 * h * i for j in range(1, MAX - 1)] for i in range(1,  MAX - 1)]
# Y = [[(j - i/2.0) * h for j in range(1, MAX - 1)] for i in range(1, MAX - 1)]
# Z = [[0.0 for j in range(1, MAX - 1)] for i in range(1, MAX - 1)]

f_array = [[0.0 for j in range(MAX)] for i in range(MAX)]
boundary_data = [[0.0 for j in range(MAX)] for i in range(MAX)]

ORDER = [[-1 for j in range(MAX)] for i in range(MAX)]
MAP = []


def initialize_all():
	global ORDER, MAP, N, MAX, g, f_array, boundary_data, S3
	IND = 0

	for i in range(1, N + 1):
		for j in range(1, N + i):
			ORDER[i][j] = IND
			MAP.append((i, j))
			IND += 1

	for i in range(N + 1, 2*N):
		for j in range(i - N + 1, 2*N):
			ORDER[i][j] = IND
			MAP.append((i, j))
			IND += 1

	g = [0.0 for k in range(IND)]

	for i in range(MAX):
		for j in range(MAX):
			x = S3 * h * i
			y = (j - i/2.0) * h

			# f_array[i][j] = det(x, y, h)
			# f_array[i][j] = F(x, y, h)
			# f_array[i][j] = math.cos(x) * math.cos(y) + 1
			# f_array[i][j] = 10
			# if i > N:
			# 	f_array[i][j] = 0.1
			if i < N/2 or i > 3*N/2:
				f_array[i][j] = 1
			else:
				f_array[i][j] = 0.1
			# f_array[i][j] = 1
			boundary_data[i][j] = 0.0
			# boundary_data[i][j] = F(x - 10, y - 5)



def val(vg, i, j):
	global ORDER, boundary_data
	ind = ORDER[i][j]
	if ind == -1:
		return boundary_data[i][j]
	else:
		return vg[ind]


def next_step(vg):
	global ORDER, MAP, f_array
	N = len(vg)

	res = []
	for k in range(N):
		(i, j) = MAP[k]

		a_d = val(vg, i + 1, j) + val(vg, i - 1, j)
		b_e = val(vg, i + 1, j + 1) + val(vg, i - 1, j - 1)
		c_f = val(vg, i, j + 1) + val(vg, i, j - 1)

		A = 4
		B = -4*(a_d + b_e + c_f) / 3.0
		C = ((a_d + b_e + c_f)**2  - 2*(a_d**2 + b_e**2 + c_f**2))/3.0 - f_array[i][j]

		if B**2 - 4*A*C <= 0:
			res.append(-B/(2*A))
		else:
			res.append((-B - math.sqrt(B**2 - 4*A*C))/(2*A))

		res[k] = (res[k] + 9*vg[k])/10.0

	return res


def update(event):
	global X, Y, Z, g
	for i in range(20):
		g = next_step(g)

	# Z = np.array([[0.0 for j in range(1, MAX - 1)] for i in range(1, MAX - 1)])
	# for i in range(1, MAX - 1):
	# 	for j in range(1, MAX - 1):
	# 		Z[i - 1][j - 1] += (val(g, i, j) - val(g, i + 1, j))**2 / 6.0
	# 		Z[i - 1][j - 1] += (val(g, i, j) - val(g, i + 1, j + 1))**2 / 6.0
	# 		Z[i - 1][j - 1] += (val(g, i, j) - val(g, i, j + 1))**2 / 6.0
	# 		Z[i - 1][j - 1] += (val(g, i, j) - val(g, i - 1, j))**2 / 6.0
	# 		Z[i - 1][j - 1] += (val(g, i, j) - val(g, i - 1, j - 1))**2 / 6.0
	# 		Z[i - 1][j - 1] += (val(g, i, j) - val(g, i, j - 1))**2 / 6.0

	Z = np.array([[val(g, i, j) for j in range(MAX)] for i in range(MAX)])
	ax.clear()
	ax.plot_surface(X, Y, Z, cmap=cm.coolwarm, linewidth=0, antialiased=False)
	plt.draw()

# ========================================================================================================================


fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
ax.autoscale(enable=None, axis="x", tight=True)
ax.set_ylim(-(MAX - 1)*h/4.0, (MAX - 1)*h*3.0/4.0)
ax.set_xlim((MAX - 1)*h*(S3 - 1)/2, (MAX - 1)*h*(S3 + 1)/2)

initialize_all()
Z = np.array([[val(g, i, j) for j in range(MAX)] for i in range(MAX)])
# Z = np.array([[val(g, i, j) for j in range(1, MAX - 1)] for i in range(1, MAX - 1)])
# Z = np.array([[f_array[i][j] for j in range(1, MAX - 1)] for i in range(1, MAX - 1)])

surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm, linewidth=0, antialiased=False)

axbutton = fig.add_axes([0.81, 0.05, 0.1, 0.075])
button = Button(axbutton, 'iter')
button.on_clicked(update)


plt.show()

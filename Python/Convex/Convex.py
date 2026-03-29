# A program that tries to iteratively solve the Monge-Ampere equation.

import math
import scipy
import random

def F(x, y):
	return (2*x**2 + 3*y**2 - 4*x*y)/2.0 + 0.529729 + 0.9472834982 * x - 0.354235239 * y

def det(x, y, h):
	S = math.sqrt(3)/2.0
	w = F(x, y)
	a = F(x, y + h)
	b = F(x + S*h, y + h/2)
	c = F(x + S*h, y - h/2)
	d = F(x, y - h)
	e = F(x - S*h, y - h/2)
	f = F(x - S*h, y + h/2)

	return ((a + b + c + d + e + f - 2*w)**2 + 8*w**2 - 2*((a + d)**2 + (b + e)**2 + (c + f)**2))/(3.0 * h**4)
	# return (2*(a + d + b + e + c + f) - 12*w)/3.0

# for i in range(10):
# 	h = 2**(-i)
# 	print(det(0, 0, 2**(-i)) * 4**i)
# 	# print(det(0, 0, 2**(-i)) * 16**i)
# 	# print(math.log2(0.00000000000001 + abs(det(0, 0, 2**(-i)))))

# So error is of order O(h^2)



N = 2
MAX = 100
h = 1
ORDER = [[-1 for i in range(MAX)] for j in range(MAX)]
# f_array = [[det(math.sqrt(3)/2 * h * i, (j - i/2.0) * h, h) for i in range(MAX)] for j in range(MAX)]
# boundary_data = [[F(math.sqrt(3)/2 * h * i, (j - i/2.0) * h) for i in range(MAX)] for j in range(MAX)]
f_array = [[1 for i in range(MAX)] for j in range(MAX)]
boundary_data = [[0 for i in range(MAX)] for j in range(MAX)]
# (i, j) corresponds to i * (sqrt(3)/2, -1/2) + j * (0, 1)
MAP = []

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

# Now ORDER[MAP[k]] = k

# Next, initialize f_array

for k in range(IND):
	(i, j) = MAP[k]

def val(vg, i, j):
	global ORDER, boundary_data
	ind = ORDER[i][j]
	if ind == -1:
		return boundary_data[i][j]
	else:
		return vg[ind]

def func(vg):
	global ORDER, MAP, f_array

	res = []
	last = 0.0
	for k in range(len(vg) - 1):
		(i, j) = MAP[k]
		w = vg[k]

		a_d = val(vg, i + 1, j) + val(vg, i - 1, j)
		b_e = val(vg, i + 1, j + 1) + val(vg, i - 1, j - 1)
		c_f = val(vg, i, j + 1) + val(vg, i, j - 1)

		res.append(((a_d + b_e + c_f - 2*w)**2 + 8*w**2 - 2*(a_d**2 + b_e**2 + c_f**2))/3.0 - f_array[i][j])

		last += min(0.0, (2*(a_d + b_e + c_f) - 12*w)/3.0)
	
	res.append(last**2 + vg[-1]**2)
	# print(vg)
	# print(res)
	return res

def jac(vg):
	global ORDER, MAP
	N = len(vg)

	res = [[0.0 for l in range(N)] for k in range(N)]
	last = 0.0
	for k in range(N - 1):
		(i, j) = MAP[k]

		w = vg[k]
		a_d = val(vg, i + 1, j) + val(vg, i - 1, j)
		b_e = val(vg, i + 1, j + 1) + val(vg, i - 1, j - 1)
		c_f = val(vg, i, j + 1) + val(vg, i, j - 1)

		# ((a_d + b_e + c_f - 2*w)**2 + 8*w**2 - 2*(a_d**2 + b_e**2 + c_f**2))/3.0
		tr = (2*(a_d + b_e + c_f) - 12*w)/3.0
		dtr = 0.0
		if tr < 0.0:
			dtr = 2.0/3.0
			res[N - 1][k] += -4
			last += tr

		S = 2*(a_d + b_e + c_f - 2*w)
		res[k][k] = (-2*S + 16*w)/3.0

		l = ORDER[i + 1][j]
		if l != -1:
			res[k][l] = (S - 4*a_d)/3.0
			res[N - 1][l] += dtr
		l = ORDER[i - 1][j]
		if l != -1:
			res[k][l] = (S - 4*a_d)/3.0
			res[N - 1][l] += dtr
		l = ORDER[i + 1][j + 1]

		if l != -1:
			res[k][l] = (S - 4*b_e)/3.0
			res[N - 1][l] += dtr
		l = ORDER[i - 1][j - 1]
		if l != -1:
			res[k][l] = (S - 4*b_e)/3.0
			res[N - 1][l] += dtr

		l = ORDER[i][j + 1]
		if l != -1:
			res[k][l] = (S - 4*c_f)/3.0
			res[N - 1][l] += dtr
		l = ORDER[i][j - 1]
		if l != -1:
			res[k][l] = (S - 4*c_f)/3.0
			res[N - 1][l] += dtr

	for i in range(N - 1):
		res[N - 1][i] *= 2*last
	res[N - 1][N - 1] = 2*vg[N - 1]

	return res

# print(jac([random.random() for _ in range(IND)]))


# print(jac([0, 0, 0, 0, 0, 0, 0, 0]))


print(scipy.optimize.root(func, [-1.0 for _ in range(IND + 1)], jac=jac))
# print(scipy.optimize.root(func, [-1.22076, -1.22076, -1.22076, -1.72076, -1.22076, -1.22076, -1.22076, 0], jac=jac))
# print(scipy.optimize.root(func, [-1.32076, -1.12076, -1.22076, -1.92076, -1.02076, -1.22076, -1.32076, 0], jac=jac))

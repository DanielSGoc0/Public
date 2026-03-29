# Basically, looks for an optimal lattice on which to perform integration
# This program analyzes the integration schemes in probabilistic setting.

from tkinter import *
import tkinter
import copy
import random
import math
import cmath
import numpy as np
import scipy
from itertools import product

# W = 19
# H = 22

# L = 21
D = 1
N = 6
L = N
R = 10 # radius, the size up to which we consider the elements

X = [] # These are points

def evaluate_sum():
	global L, D, R
	T = R // L
	SUM = -1.0

	for a in product(range(-T, T + 1), repeat=D):
		M = 0.0
		for v in a:
			M = M + v*v
		SUM = SUM + math.exp(-2*L*math.sqrt(M))

	return SUM

# 	return SUM/math.sqrt(2)

# def evaluate_sum2():
# 	global L, D, R, W, H
# 	SUM = -1.0

# 	for i in range(-R, R + 1, 1):
# 		for j in range(-R, R + 1, 1):
# 			if j % W == 0:
# 				if (2*i + H*(j/W)) % (2*H) == 0:
# 					SUM = SUM + math.exp(-math.sqrt(i*i + j*j))

# 	return SUM/math.sqrt(2)

def generate_random_points():
	global D, N, R, X
	X = []

	for _ in range(N):
		P = []
		for _ in range(D):
			P.append(random.random() * 2*math.pi)
		X.append(P)

def generate_normal_points():
	global D, N, L, X
	X = []

	for a in product(range(L), repeat=D):
		vec = []
		for i in range(D):
			vec.append((2*a[i] + 1) * 2*math.pi / (2*L))
		X.append(vec)

# def generate_triangular_points():
# 	global X, W, H
# 	for i in range(H):
# 		for j in range(W):
# 			X.append(((i*2*math.pi)/H, (2*math.pi*(2*j + i))/(2*W)))
	

def calculate_problem():
	global D, N, R, X
	
	A = []
	e = []

	for a in product(range(-R, R + 1), repeat=D):
		M = 0
		for v in a:
			M = M + v*v

		if R*R >= M:
			if a[0] % 2 == 1:
				e.append(1/a[0])
			else:
				e.append(0)
			M = math.sqrt(M)

			vec = []
			for i in range(N):
				S = 0
				for j in range(D):
					S = S + a[j] * X[i][j]
				vec.append(cmath.exp(complex(-M, S)))
			
			A.append(vec)

	A = np.transpose(np.array(A))
	# print(A)

	D1 = np.matmul(A, np.conjugate(np.transpose(A)))
	# print(D1)
	# print(np.linalg.det(D1))

	A = np.append(A, [e], axis=0)
	# print(A)

	D2 = np.matmul(A, np.conjugate(np.transpose(A)))
	# print(D2)
	# print(np.linalg.det(D2))

	res = abs(np.linalg.det(D2))/abs(np.linalg.det(D1))
	print(res)
	return res

	# print(res.x)
	# print(res.con)

# generate_random_points()
# generate_triangular_points()
# print(evaluate_sum())
# print(evaluate_sum2())
# calculate_problem()
# print(optimize_problem())

OPT = 1000
for i in range(3000):
	print("============  i = " + str(i))
	generate_random_points()
	# generate_normal_points()
	OPT = min(OPT, calculate_problem())

print("====================================")
generate_normal_points()
# print(evaluate_sum())
calculate_problem()
# print(math.exp(-2*N))
print(OPT)
	
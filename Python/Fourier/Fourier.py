# Basically, looks for an optimal lattice on which to perform integration

from tkinter import *
import tkinter
import copy
import random
import math
import numpy
import scipy
from itertools import product

W = 19
H = 22

L = 21
D = 2
N = W*H
R = 32 # radius, the size up to which we consider the elements

X = [] # These are points

def evaluate_sum():
	global L, D, R
	T = R // L
	SUM = -1.0

	for a in product(range(-T, T + 1), repeat=D):
		M = 0.0
		for v in a:
			M = M + v*v
		SUM = SUM + math.exp(-L*math.sqrt(M))

	return SUM/math.sqrt(2)

def evaluate_sum2():
	global L, D, R, W, H
	SUM = -1.0

	for i in range(-R, R + 1, 1):
		for j in range(-R, R + 1, 1):
			if j % W == 0:
				if (2*i + H*(j/W)) % (2*H) == 0:
					SUM = SUM + math.exp(-math.sqrt(i*i + j*j))

	return SUM/math.sqrt(2)

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
			vec.append(a[i] * 2*math.pi / L)
		X.append(vec)

def generate_traingular_points():
	global X, W, H
	for i in range(H):
		for j in range(W):
			X.append(((i*2*math.pi)/H, (2*math.pi*(2*j + i))/(2*W)))
	

def optimize_problem():
	global D, N, R, L, X
	
	A_eq = []
	b_eq = []
	bounds = []
	c = []
	x0 = []

	c.append(1)
	x0.append(0)
	A_eq.append(numpy.ones(N).tolist())
	bounds.append((None, None))

	for a in product(range(-R, R + 1), repeat=D):
		value = 0
		for v in a:
			if v > 0:
				value = 1
				break
			elif v < 0:
				value = -1
				break

		if value == 1:
			M = 0
			for v in a:
				M = M + v*v
			M = math.sqrt(M)
			M = math.sqrt(2) * math.exp(-M)

			x0.append(0)
			x0.append(0)
			c.append(0)
			c.append(0)
			bounds.append((-M, M))
			bounds.append((-M, M))

			vec1 = []
			vec2 = []
			for i in range(N):
				S = 0
				for j in range(D):
					S = S + a[j] * X[i][j]
				vec1.append(math.cos(S))
				vec2.append(math.sin(S))
			
			A_eq.append(vec1)
			A_eq.append(vec2)

	for i in range(N):
		b_eq.append(0)

	A_eq = numpy.array(A_eq)
	A_eq = numpy.transpose(A_eq)
	# print(A_eq)

	res = scipy.optimize.linprog(c, x0 = x0, A_eq = A_eq, b_eq = b_eq, bounds = bounds)

	print(evaluate_sum())
	print(-res.fun)
	if -res.fun < evaluate_sum() * (1 - 0.001):
		print("SHIIIIIIIIIIITTT!!!!")
	return -res.fun
	# print(res.x)
	# print(res.con)

generate_traingular_points()
print(evaluate_sum())
print(evaluate_sum2())
optimize_problem()
# print(optimize_problem())

# OPT = 1000
# for i in range(3000):
# 	print("============  i = " + str(i))
# 	generate_random_points()
# 	# generate_normal_points()
# 	OPT = min(OPT, optimize_problem())
# print("====================================")
# print(evaluate_sum())
# print(OPT)
	
# A generic old inequality checker

from tkinter import *
import tkinter
import copy
import random
import math
import numpy
import scipy

a = 1.1

def f(x):
	global a
	# return math.log(x)
	# return -x*math.log(x)
	# return -x*math.log(x) - (1 - x)*math.log(1 - x)
	# return x**a + (1 - x)**a
	# return x**a
	# return -x**a
	return - x**a - (1 - x)**a
	# if x < 1.0/3:
	# 	return (x - 1.0/3) * 10
	# elif x < 2.0/3:
	# 	return 0
	# else:
	# 	return -(x - 2.0/3) * 1

def df(x):
	global a
	# return 1/x
	# return -math.log(x) - 1
	# return math.log((1 - x)/x)
	# return a*x**(a - 1) - a*(1-x)**(a - 1)
	# return a*x**(a - 1)
	# return -a*x**(a - 1)
	return -a * x**(a - 1) + a *(1-x)**(a - 1)
	# if x < 1.0/3:
	# 	return 10
	# elif x < 2.0/3:
	# 	return 0
	# else:
	# 	return -1

def ddF(x):
	global a
	# return -1/(x*x)
	# return -1/x
	# return -1/(x*(1-x))
	# return a*(a - 1)*x**(a - 2) + a*(a - 1)*(1-x)**(a - 2)
	# return a*(a - 1)*x**(a - 2)
	# return -a*(a - 1)*x**(a - 2)
	return -a*(a - 1)*x**(a - 2) - a*(a - 1)*(1-x)**(a - 2)

def optimal_m(v):
	a = 0.0
	b = 1.0
	for _ in range(70):
		h = (a + b)/2
		if df(h) > v:
			a = h
		else:
			b = h
	return (a + b)/2


def optimum(a, b):
	S = (f(b) - f(a))/(b - a)
	m = optimal_m(S)
	# print(m - 1/(1 + math.exp(S)))
	return f(m) - ((b - m) * f(a) + (m - a) * f(b))/(b - a)
	# return m

def optimal_vals(a, b, c):
	df1 = (f(b) - f(a))/(b - a)
	df2 = (f(c) - f(b))/(c - b)
	df3 = (f(c) - f(a))/(c - a)
	dfa = df(a)
	dfb = df(b)
	dfc = df(c)
	A = [[0 for i in range(3)] for j in range(3)]
	B = [0 for _ in range(3)]
	A[0][0] = -(dfa - df1)/(b - a)
	A[0][1] = 0
	A[0][2] = (dfa - df3)/(c - a)
	B[0] = - b*(dfa - df1)/(b - a) + c*(dfa - df3)/(c - a)
	A[1][0] = (dfb - df1)/(b - a)
	A[1][1] = -(dfb - df2)/(c - b)
	A[1][2] = 0
	B[1] = a*(dfb - df1)/(b - a) - c*(dfb - df2)/(c - b)
	A[2][0] = 0
	A[2][1] = (dfc - df2)/(c - b)
	A[2][2] = -(dfc - df3)/(c - a)
	B[2] = b*(dfc - df2)/(c - b) - a*(dfc - df3)/(c - a)
	I = scipy.linalg.solve(A, B)
	return I

def check_hypothesis(a, b, c):
	I = optimal_vals(a, b, c)
	if a < I[0] and I[0] < b and b < I[1] and I[1] < c and I[0] < I[2] and I[2] < I[1]:
		return True
	else:
		return False

def difference(a, b, c):
	# print(optimum(a, c))
	# print(optimum(a, b))
	# return math.sqrt(optimum(a, c)) - (math.sqrt(optimum(a, b)) + math.sqrt(optimum(b, c)))
	return optimum(a, c) - (optimum(a, b) + optimum(b, c))

def second_derivative(a, b):
	I = optimal_m((f(b) - f(a))/(b - a))
	R = (df(b) - df(I))/(b - a)
	T = (I - a)/(b - a)

	return -(R*R/ddF(I) - 2*R*T + T*ddF(b))
	# return R*R/ddF(I) + T*ddF(b)

# print(optimal_vals(0.000001, 1.0/2, 0.999999))

# for i in range(1000000):
# 	if i % 1000 == 0:
# 		print(i)
# 	tab = [random.random() for _ in range(3)]
# 	tab = sorted(tab)
# 	# print(tab)
# 	# if(check_hypothesis(tab[0], tab[1], tab[2])):
# 		# print("OK")
# 		# pass
# 	# else:
# 		# print("NOT OK")
# 	res = difference(tab[0], tab[1], tab[2])
# 	if res < -0.0001:
# 		print(res)
# 		print(tab)
# 	# else: 
# 	# 	print("OK")

def I(a, b):
	return 1/(1 + math.exp((f(b) - f(a))/(b - a)))

def distance(a, I, b):
	return f(I) - ((I - a) * f(b) + (b - I) * f(a))/(b - a)

def weighted(a, I1, b, I2, c):
	# gamma = (b - a)/(c - a)
	# I3 = gamma*I1 + (1 - gamma)*I2

	return second_derivative(I1, I2)
	# return difference(a, b, c)

	# if I3 < 0.000001 or I3 > 0.9999999:
	# 	return 0
	# PART1 = f(I1) - ((I1 - a) * f(b) + (b - I1) * f(a))/(b - a)
	# PART2 = f(I2) - ((I2 - b) * f(c) + (c - I2) * f(b))/(c - b)
	# PART3 = f(I3) - ((I3 - a) * f(c) + (c - I3) * f(a))/(c - a)
	# PART1 = f(I1) - (I1 * f(b) - I1 * f(a))/(b - a)
	# PART2 = f(I2) - (I2 * f(c) - I2 * f(b))/(c - b)
	# PART3 = f(I3) - (I3 * f(c) - I3 * f(a))/(c - a) + b
	# return PART3 - (PART1 + PART2)

	# X1 = f(I1) + (a - I1) * (f(b) - f(a))/(b - a)
	# X2 = f(I2) + (c - I2) * (f(c) - f(b))/(c - b)
	# X3 = f(I3) + (b - I3) * (f(c) - f(a))/(c - a) + gamma * f(a) + (1 - gamma) * f(c)
	# return (X3 - (X1 + X2))

	# return a*b*math.log(a/b)/(a - b) + b*c*math.log(b/c)/(b - c) - a*c*math.log(a/c)/(a - c) - b

	# print(optimal_m((f(b) - f(a))/(b - a)) - I(a, b))
	# return 0
	# return (optimum(b, c) - b)/(c - b) - (optimum(a, b) - a)/(b - a)
	# return (optimum(b, c) - b)/(c - b) - (optimum(a, c) - a)/(c - a)
	# return (optimum(a, c) - a)/(c - a) - (optimum(a, b) - a)/(b - a)
	# return  I(b, c) - b - I(a, b) - I(a, c)
	# S = c * math.log(1/(c * I2) - 1)
	# T = b * math.log(1/(b * I1) - 1)
	# if S <= math.log(2) and T <= math.log(2) and 1/math.e <= I1 and I2 <= 1/2 and T/b - S/c >= 2.45*(c - b) and S/(1 - c) - T/(1 - b) >= 2.45*(c - b) and T >= 4*b*(1 - b) and S >= 4*c*(1-c) and (S - T)/(c - b) >= math.log((1 - c)/c) and (S - T)/(c - b) <= math.log((1 - b)/b):
	# 	return -(S - T - (c - b) * math.log(1/(b + (c - b)*I1) - 1))
	# else:
	# 	return 0
	# d = (f(b) - f(a))/(b - a)
	# e = math.exp(d)
	# w = (-b*math.log(a/b) - (1 - b)*math.log((1-a)/(1-b)))/(b-a)**2
	# print(w * e/(1+e)**2)
	# return w * e/(1+e)**2 - (b - 1/(1 + e))/(b - a)
	# return 1
	# A = 1/(b * (1 + math.exp(f(b)/b)))
	# A = (1 + random.random())/2.0
	# B = b
	# x = c - b

	# return (A*x)/((A*x + B - 1)*(A*x + B)) - math.log((A*x + B)/(x + B)) - ((1 - B - x)/(1 - A*x - B) - 1)
	# return -A/((1 - A*x - B)*(A*x + B)) + (1 - A)/(x + B) * (A*x + x + 2*B)/(A*x + B) * 1.0/2 + (1 - A)/(1 - A*x - B)
	# return -2*A*(x + B) + (1 - A)*(x + B)*(1 + A*x + B) + (1 - A)*(A*x + B)*(1 - A*x - B)
	# return B*(1 - A) - A
	# print((I(a, c) - a)/(c - a))
	# print((I(b, c) - b)/(c - b))
	# print()
	# return (I(a, b) - a)/(b - a)


# for i in range(1000000):
# 	if i % 1000 == 0:
# 		print(i)
# 	tab = [random.random() for _ in range(5)]
# 	tab = sorted(tab)
# 	# print(tab)
# 	# if(check_hypothesis(tab[0], tab[1], tab[2])):
# 		# print("OK")
# 		# pass
# 	# else:
# 		# print("NOT OK")
# 	res = weighted(tab[0], tab[1], tab[2], tab[3], tab[4])
# 	# print(res)
# 	if res < -0.0001:
# 		print(res)
# 		print(tab)
# 	# else: 
# 	# 	print("OK")

def check_interval(a, b):
	tab1 = [a + random.random() * (b - a) for _ in range(5)]
	tab1 = sorted(tab1)
	tab2 = [a + random.random() * (b - a) for _ in range(2)]
	res = weighted(tab1[0], tab2[0], tab1[2], tab2[1], tab1[4])
	if res < -0.000001:
		print(res)
		print(tab1)
		print(tab2)

def check_intervals():
	N = 1000
	h = 1.0/N
	for i in range(N):
		print(i)
		for _ in range(1000):
			check_interval(i*h, (i + 1)*h)

# check_intervals()

# for i in range(10000000):
# 	mu = random.random()
# 	gamma = random.random()

# 	tab = [random.random() for _ in range(2)]
# 	tab = sorted(tab)
# 	xa = tab[0]
# 	xc = tab[1]
# 	xb = gamma * xa + (1 - gamma) * xc	

# 	tab = [random.random() for _ in range(2)]
# 	tab = sorted(tab)
# 	ya = tab[0]
# 	yc = tab[1]
# 	yb = gamma * ya + (1 - gamma) * yc
	
# 	fa = f(xa) + f(ya) - f(mu * xa + (1 - mu) * ya)
# 	fb = f(xb) + f(yb) - f(mu * xb + (1 - mu) * yb)
# 	fc = f(xc) + f(yc) - f(mu * xc + (1 - mu) * yc)

# 	diff = fb - gamma * fa - (1 - gamma) * fc
# 	if diff < 0:
# 		print("OK")
# 	else:
# 		print("NOT OK")
# 		# pass

# for i in range(10000000):
# 	if i % 1000 == 0:
# 		print(i)
# 	tab1 = [random.random() for _ in range(5)]
# 	tab1 = sorted(tab1)
# 	# tab1[1] = tab1[0]*0.999999 + tab1[1] * (1 - 0.999999)
# 	tab2 = [random.random() for _ in range(2)]
# 	tab2 = sorted(tab2)
# 	# print(tab)
# 	# if(check_hypothesis(tab[0], tab[1], tab[2])):
# 		# print("OK")
# 		# pass
# 	# else:
# 		# print("NOT OK")
# 	# tab1[0] = 0.000001
# 	# res = weighted(tab1[0], tab1[1], tab1[2], tab1[3], tab1[4])
# 	res = weighted(tab1[0], tab2[0], tab1[2], tab2[1], tab1[4])
# 	# print(res)
# 	# print(res)
# 	if res < -0.000001:
# 		print(res)
# 		print(tab1)
# 		print(tab2)
# 	# else: 
# 	# 	print("OK")

A = [[random.random() for j in range(2)] for i in range(2)]
A = numpy.array(A)

B = [[random.random() for j in range(2)] for i in range(2)]
B = numpy.array(B)

print(numpy.linalg.det(A + B) + numpy.linalg.det(A - B))
print(2*numpy.linalg.det(A) + 2*numpy.linalg.det(B))

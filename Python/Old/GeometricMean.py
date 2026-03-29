# Honestly, I don't know what this one is about

import random
import math
import numpy as np
from scipy import linalg

def GeometricMean2(a, b, c, d, e, f):
	A = np.array([[a, b], [b, c]])
	B = np.array([[d, e], [e, f]])

	Asqrt = linalg.sqrtm(A)
	Asqrtinv = linalg.inv(Asqrt)
	W = np.dot(Asqrtinv, np.dot(B, Asqrtinv))
	W = linalg.sqrtm(W)
	W = np.dot(Asqrt, np.dot(W, Asqrt))

	return [W[0][0], W[0][1], W[1][1]]


def GeometricMean(a, b, c, d, e, f):

	det1 = a*c - b**2
	det2 = d*f - e**2

	A = a*e - b*d
	B = a*f - c*d
	C = b*f - c*e

	DET = B**2 - 4*A*C

	x = (-B + math.sqrt(DET))/(2*A)
	z = (-B - math.sqrt(DET))/(2*A)

	d1 = math.sqrt((a*x**2 + 2*b*x + c)*(d*x**2 + 2*e*x + f))
	d2 = math.sqrt((a*z**2 + 2*b*z + c)*(d*z**2 + 2*e*z + f))

	return A**2/DET * np.array((d1 + d2, -d1*z - d2*x, d1*z**2 + d2*x**2))

def random_positive():
	while True:
		a = random.random()
		b = random.random() * 2 - 1
		c = random.random()

		if a*c - b**2 > 0 and a + c > 0:
			return (a, b, c)


# (a, b, c) = random_positive()
# (d, e, f) = random_positive()

# (a, b, c) = (2, 0, 4)
# (d, e, f) = (1, -1, 3)
# (a, b, c) = (2, 1, 2)
# (d, e, f) = (3, 2, 3)
# print((a, b, c))
# print((d, e, f))

# print(GeometricMean(a, b, c, d, e, f))
# print(GeometricMean2(a, b, c, d, e, f))


def Hessian1(x, y):
	# f(x, y) = x^2 - x*y + y^2

	return (2, -1, 2)


def Hessian2(x, y):
	# g(x, y) = e^(3x - 2y)
	v = math.exp(3*x - 2*y)
	w = math.exp(x + y)

	return (w + 9*v, w - 6*v, w + 4*v)

def EvalMean(x, y):
	(a, b, c) = Hessian1(x, y)
	(d, e, f) = Hessian2(x, y)
	return GeometricMean(a, b, c, d, e, f)

for k in range(30):
	h = 2**(-k)

	x = 0.8
	y = 1.4
	
	(a1, b1, c1) = EvalMean(x, y)
	(a2, b2, c2) = EvalMean(x, y + h)

	print((a2 - a1)/h, (b2 - b1)/h, (c2 - c1)/h)

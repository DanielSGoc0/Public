# NOT concave.
import math
import random

def f(x):
	# return x*(1-x)
	# return math.cos(x)
	# return math.sqrt(x)
	# return -1.0/(x + 1.0)
	return -math.tan(x)

def S(x, n):
	r = 0.0
	for i in range(0, n):
		r += (f(x[i + 1]) + f(x[i])) * (x[i + 1] - x[i])
	return r

n = 4
p = 0.5

for attempt in range(100000):
	if attempt % 1000 == 0:
		print(attempt)

	x = [random.random() for i in range(n - 1)]
	x += [0.0, 1.0]
	x = sorted(x)

	y = [random.random() for i in range(n - 1)]
	y += [0.0, 1.0]
	y = sorted(y)

	z = [p*x[i] + (1.0-p)*y[i] for i in range(n + 1)]

	w = S(z, n) - p*S(x, n) - (1.0 - p)*S(y, n)
	if w < 0.0:
		print(x)
		print(y)
		print(z)
		print(w)
		break

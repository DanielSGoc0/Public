# A generic tester of specific inequality

import random

def val(x, y, z):
	return (1.0/(x + y + z) + 1.0/(1.0/x + 1.0/y + z) + 1.0/(1.0/x + y +1.0/z) + 1.0/(x + 1.0/y + 1.0/z))

def val2(x, y, z):
	return val(x, y, z) * val(1.0/x, 1.0/y, 1.0/z)

def LHS(x, y, z):
	return (x + y + z)*(x**2 + y**2 + z**2)

def RHS(x, y, z):
	return 1.0/15.0 * (1.0/x + 1.0/y + 1.0/z) * ((x + y)**4 + (y + z)**4 + (z + x)**4 - x**4 - y**4 - z**4)

MIN = 10000000000
MAX = 0
for k in range(1000000):
	a = 1.0/random.random() - 1.0
	b = 1.0/random.random() - 1.0
	c = 1.0/random.random() - 1.0

	# v = val2(a, b, c)
	# MIN = min(MIN, v)
	# MAX = max(MAX, v)
	print(LHS(a, b, c), RHS(a, b, c))

print(MIN)
print(MAX)

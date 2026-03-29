# Used to check certain numeric-theoretical hypothesis about representability
# of two sequences as products of Binomial Coefficients.

from random import *

def GCD(x, y):
	if x == 0:
		return y
	if y == 0:
		return x
	if x < y:
		return GCD(x, y%x)
	return GCD(y, x%y)


# check if two sequences are "cool2". We assume sum_i a_i = sum_j b_j
def checker2(a, b):
	
	for v in a:
		for k in range(1, v):
			# we set x = k/v
			SUM1 = 0
			SUM2 = 0

			for x in a:
				# calculate floor(x * k/v - epsilon)
				SUM1 += (x*k - 1)//v

			for y in b:
				# calculate floor(y * k/v - epsilon)
				SUM2 += (y*k - 1)//v
			
			if SUM2 > SUM1:
				return False

	return True
		

# check if two sequences are "cool1". We assume sum_i a_i = sum_j b_j
def checker1(a, b):
	if len(a) == 0:
		return True
	
	for i in range(len(a)):
		if b[0] == a[i]:
			c = []
			for i2 in range(i):
				c.append(a[i2])
			for i2 in range(i + 1, len(a)):
				c.append(a[i2])

			d = []
			for j in range(1, len(b)):
				d.append(b[j])

			return checker1(c, d)

		elif b[0] < a[i]:
			c = []
			for i2 in range(i):
				c.append(a[i2])
			c.append(a[i] - b[0])
			for i2 in range(i + 1, len(a)):
				c.append(a[i2])

			d = []
			for j in range(1, len(b)):
				d.append(b[j])

			if checker1(c, d):
				return True
	return False
		

for i in range(1000):
	SUM = randint(10, 20)

	a = []
	SUM1 = 0
	while SUM1 < SUM:
		x = randint(1, 7)
		if SUM1 + x > SUM:
			a.append(SUM - SUM1)
			break

		a.append(x)
		SUM1 += x

		if SUM1 == SUM:
			break

	b = []
	SUM2 = 0
	while SUM2 < SUM:
		y = randint(1, 4)
		if SUM2 + y > SUM:
			b.append(SUM - SUM2)
			break

		b.append(y)
		SUM2 += y

		if SUM2 == SUM:
			break

	cool1 = checker1(a, b)
	cool2 = checker2(a, b)

	print(cool1, cool2)
	if cool1 != cool2:
		print(a)
		print(b)
		break

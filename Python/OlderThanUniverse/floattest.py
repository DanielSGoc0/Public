
from bigfloat import *


def Binomial(a, b):
	if b > a:
		return 0
	if b == 0:
		return 1
	if b > a/2:
		return Binomial(a, a - b)
	w = 1
	for i in range(1, b + 1):
		w *= a + 1 - i
		w /= i
	return w

def funkcja1(x):
	N = 13
	if x % 2 == 0:
		return x
	else:
		# MAIN PART
		# argument to obecne x
		with precision(200):
			y = BigFloat(0)
			#sig = BigFloat(1)
			for k in range(0, N + 1):
				z = BigFloat(k + 1)
				z = mul(z, mul(BigFloat(Binomial(N + k + 1, N)), BigFloat(Binomial(N, k))))
				z = mul(z, div(BigFloat(N + 1 - x), BigFloat(k + 1 - x)))
				#print z
				#z = mul(sig, z)
				#y = add(y, z)
				#sig = neg(sig)
				if (N - k) % 2 == 0:
					y = add(y, z)
				else:
					y = sub(y, z)
				#print y
			#print y
			for i in range(1, N + 1):
				y = mul(y, div(BigFloat(x - i), BigFloat(x + i)))
			# MAIN PART
			return y

def funkcja2(x):
	N = 13
	if x % 2 == 0:
		return x
	else:
		q = x
		with precision(200):
			# MAIN PART
			# argument to obecne x
			y = BigFloat(0)
			for k in range(0, N + 1):
				z = BigFloat(k + 1)
				z = mul(z, mul(BigFloat(Binomial(N + k + 1, N)), BigFloat(Binomial(N, k))))
				z = mul(z, div(BigFloat(N + 1 - x), BigFloat(k + 1 - x)))
				
				if (N - k) % 2 == 0:
					y = add(y, z)
				else:
					y = sub(y, z)
			for i in range(1, N + 1):
				y = mul(y, div(BigFloat(x - i), BigFloat(x + i)))
			# MAIN PART
			if q < 202 and q > 200:
				print y
			return y

print funkcja1(15)
import math
import scipy

for p in range(1, 30):
	c = 0.0
	w = 0.0

	for k in range(1, p + 1):
		c += (-1)**k * scipy.special.binom(p, k) * math.log(k)
		w += (-1)**k * scipy.special.binom(p, k) * math.log(k)**2

	res = w + c**2

	print(res)

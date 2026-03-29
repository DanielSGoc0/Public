# Suppose that ad | b^2 - c^2 and bc | a^2 - d^2
# and that a,b,c,d are all pairwise coprime.
# Must (b^2 - c^2)(a^2 - d^2) / (abcd) be a square?
import math
import random 

if False:
	for attempt in range(1000):
		a = 2 * random.randint(0, 100) + 1

		d = None
		while True:
			d = 2 * random.randint(0, 100) + 1
			if math.gcd(a, d) == 1 and (not a == d):
				break
			
		x = abs(d**2 - a**2)

		b = None
		while True:
			b = random.randint(1, x)
			if x % b == 0:
				break
		
		for c in range(1, x):
			if x % (b*c) == 0 and math.gcd(b, c) == 1 and (c**2 - b**2) % (a*d) == 0:
				s = abs((c**2 - b**2) * (a**2 - d**2)) / (a*b*c*d)
				print()
				print(attempt)
				print(s, s**0.5)
				print(a, b, c, d)

else:
	for attempt in range(1000000):
		a = random.randint(1, 1000)

		d = None
		while True:
			d = random.randint(1, 1000)
			if math.gcd(a, d) == 1 and (not a == d):
				break
			
		x = d**2 + a**2

		b = None
		while True:
			b = random.randint(1, x)
			if x % b == 0:
				break
		
		for c in range(1, x):
			if x % (b*c) == 0 and math.gcd(b, c) == 1 and (c**2 + b**2) % (a*d) == 0 and not b == c:
				s = abs((c**2 + b**2) * (a**2 + d**2)) / (a*b*c*d)
				print()
				print(attempt)
				print(s, s**0.5)
				print(a, b, c, d)



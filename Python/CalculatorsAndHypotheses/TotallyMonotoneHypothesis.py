# Checks if particular sequence is completely monotone
import math

n = 30
k = 5
r2 = 0.01
# rho2 = 0.1
# h2 = r2 / (1.0 - rho2 * math.exp(-r2))
h2 = 1.0
rho2 = (1.0 - r2) * math.exp(r2)
print("h2 =", h2)

a = [0.0 for i in range(n)]
a[0] = r2 - math.log(rho2)
a[1] = r2
for i in range(2, n):
	a[i] = h2 * (1.0 - math.exp(-a[i-1]))
for i in range(n):
	a[i] /= h2**i

print(a)

for i in range(k, n):
	sum = 0.0
	for j in range(k + 1):
		sum += (-1)**(k - j) * a[i - j] * math.comb(k, j)
	print(sum)

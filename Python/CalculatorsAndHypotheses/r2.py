# Produces iterates of certain recursive sequence.
import math

N = 1000000
r02 = 0.999
rho02 = 0.0

r2 = [0.0 for i in range(N)]
r2[0] = r02
r2[1] = r02 * (1.0 - math.exp(-r02))/(1.0 - math.exp(-r02) * rho02)

for i in range(2, N):
	r2[i] = r2[i - 1] * (1.0 - math.exp(-r2[i - 1]))/(1.0 - math.exp(-r2[i - 2]))
	if (i & (i - 1)) == 0:
		print(i, '\t', r2[i])

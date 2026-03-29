# Checks which functions [0, inf) -> [0, inf) satisfy the following:
# if xk, yk, zk >= 0 satisfy f(xk) + f(yk) >= f(zk)
# then also f(sum xk) + f(sum yk) >= f(sum zk)
import numpy as np

def f(x):
	# return np.log(np.sqrt(x) + 1)
	# return np.log(x**2 + 1)
	# return np.sqrt(x)
	# return 1.0 - np.exp(-x)
	# return x / (1.0 + x)
	# return np.log(1.0 + x)
	# return np.sqrt(x) / (1.0 + np.sqrt(x))
	# return x + np.sqrt(x)
	# return min(x, 1.0)
	return min(x, np.sqrt(x))

def f_inv(y):
	# return (np.exp(y) - 1)**2
	return np.exp(y) - 1
	# return y**2
	# return y / (1.0 - y)

INV = False
for attempt in range(1000000):
	print(attempt)
	s = np.random.random((1))
	s = s / (1.0 - s)

	if INV:
		v = np.random.random((4)) * s
		(x0, y0, x1, y1) = v
		z0 = f_inv(f(x0) + f(y0))
		z1 = f_inv(f(x1) + f(y1))
		if z0 < 0.0 or z1 < 0.0:
			continue
		if f(x0 + x1) + f(y0 + y1) < f(z0 + z1):
			print("BAD")
			print(x0, y0, z0, x1, y1, z1)
			break

	else:
		v = np.random.random((6)) * s
		(x0, y0, z0, x1, y1, z1) = v
		if f(x0) + f(y0) >= f(z0) and f(x1) + f(y1) >= f(z1) and f(x0 + x1/2) + f(y0 + y1/2) < f(z0 + z1/2):
			print("BAD")
			print(x0, y0, z0, x1, y1, z1)
			break


	# v = np.random.random((8)) * s
	# (x0, y0, z0, w0, x1, y1, z1, w1) = v
	# if f(x0) + f(y0) + f(z0) >= f(w0) and f(x1) + f(y1) + f(z1) >= f(w1) and f(x0 + x1) + f(y0 + y1) + f(z0 + z1) < f(w0 + w1):
	# 	print("BAD")
	# 	print(x0, y0, z0, w0, x1, y1, z1, w1)
	# 	break

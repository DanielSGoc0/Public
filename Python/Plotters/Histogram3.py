# A histogram

import numpy as np
import matplotlib.pyplot as plt

D = 3
N = 3
S = 1000

# ======================================================================

def calc(M, N, D):
	MAX_amt = 1
	MAX_res = 1.0
	MAX_mean = M[0]

	for l in range(1, 2**N):
		mean = np.zeros((D), dtype='double')
		bitset = [int(n) for n in bin(l)[2:].zfill(N)]

		amt = 0
		for j in range(N):
			if bitset[j] == 1:
				mean += M[j]
				amt += 1
		mean /= amt

		res = amt
		for j in range(N):
			if bitset[j] == 1:
				res -= np.dot(M[j] - mean, M[j] - mean) / D
		
		if res > MAX_res:
			MAX_res = res
			MAX_amt = amt
			MAX_mean = mean

		if l == 2**N - 1 and MAX_amt == 1:
			MAX_mean = mean


	return (MAX_res, MAX_amt, MAX_mean)

T = 10000

X = []
ZEROS = 0
DIST = 0.0
for k in range(T):
	N = 3
	D = 3

	M = np.random.normal(0.0, 1.0, (N, D))
	(MAX_res, MAX_amt, MAX_mean) = calc(M, N, D)
	DIST += np.dot(MAX_mean, MAX_mean)
	X.append(MAX_res)
	if MAX_amt == 1:
		ZEROS += 1

print(ZEROS)
print(DIST / T)


# ======================================================================

n, bins, patches = plt.hist(x=X, bins=100, color='#0504aa', alpha=0.7, rwidth=0.85)
plt.grid(axis='y', alpha=0.75)
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title('My Very Own Histogram')
plt.text(23, 45, r'$\mu=15, b=3$')
maxfreq = n.max()
# Set a clean upper y-axis limit.
# plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)

plt.show()

# ======================================================================


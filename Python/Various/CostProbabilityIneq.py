from itertools import permutations
import random
import numpy as np

N = 10
K = 3

def set_sum_size(a):
	global N, K

	covered = [0 for i in range(N)]
	for aj in a:
		covered[aj] += 1
		if aj + K < N:
			covered[aj + K] -= 1
		if aj + K > N:
			covered[0] += 1
			covered[aj + K - N] -= 1

	total = 0
	if covered[0] > 0:
		total = 1
	for i in range(1, N):
		covered[i] += covered[i - 1]
		if covered[i] > 0:
			total += 1

	return total 


def cost_single(a):
	global N, K

	c = 1.0
	for i in range(1, len(a) + 1):
		c *= set_sum_size(a[:i]) / N

	return c


def cost(a):
	global N, K

	a = np.array(a)
	if len(a) == 0:
		return 1.0

	c = 0.0
	amt = 0

	for perm in permutations(a):
		amt += 1
		c += cost_single(perm)
	
	return c / amt


def best_removal(a, x):
	global N, K

	largest = 0.0
	for j in range(len(a)):
		if (a[j] <= x and x < a[j] + K) or (a[j] <= x + N and x + N < a[j] + K):
			b = []

			for j0 in range(j):
				b.append(a[j0])
			for j0 in range(j + 1, len(a)):
				b.append(a[j0])

			largest = max(largest, cost(b))
	
	return largest


def best_cost(a):
	global N, K

	total = 0.0
	for x in range(N):
		total += best_removal(a, x)
	
	return total / N


# =========================================================================================

# We take the definition sup_{K <= |A| < N} |amt of arcs in A|/|total amt of arcs| - |A|/N
def S(a):
	global N, K

	supremum = -1.0

	# consider intervals [i, j) so ak must satisfy i <= ak <= j - K 
	for i in range(N):
		for j in range(i + K, i + N + 1):
			amt_of_arcs = 0

			for ak in a:
				if (i <= ak and ak <= j - K) or (i <= ak + N and ak + N <= j - K):
					amt_of_arcs += 1

			supremum = max(supremum, amt_of_arcs / len(a) - (j - i) / N)
	
	return supremum


def best_removal_S(a, x):
	global N, K

	smallest = 0.0
	for j in range(len(a)):
		if (a[j] <= x and x < a[j] + K) or (a[j] <= x + N and x + N < a[j] + K):
			b = []

			for j0 in range(j):
				b.append(a[j0])
			for j0 in range(j + 1, len(a)):
				b.append(a[j0])

			smallest = min(smallest, S(b))
	
	return smallest


def best_S(a):
	global N, K

	total = 0.0
	for x in range(N):
		total += best_removal_S(a, x)
	
	return total / N



N = 20
K = 7

worst_case = 10000.0
for attempt in range(100):
	print(attempt)
	a = []

	for j in range(14):
		a.append(random.randint(0, N-1))
	a = sorted(set(a))
	if set_sum_size(a) < N:
		continue
	

	S_a = S(a)
	best_S_a = best_S(a)


	# print(S_a, best_S_a)
	improvement = S_a - best_S_a
	worst_case = min(worst_case, improvement)
	# print("improvement:", S_a - best_S_a)

	# if best_a_cost < a_cost:
	# 	print("shit")
	# 	print(a)
	# 	print(best_a_cost, a_cost)
	# 	break
print(worst_case)

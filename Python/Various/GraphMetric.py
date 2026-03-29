import numpy
import random

N = 7
powers = [2**i for i in range(N)]

def generate_random():
	global N, powers
	
	G = [[False for j in range(N)] for i in range(N)]
	for i in range(N):
		for j in range(i):
			# if random.random() < 0.5:
				G[i][j] = True
				G[j][i] = True
	return G

def calculate_metric(G):
	global N, powers

	metric = -1000000.0
	total_edges = 0
	for i in range(N):
		for j in range(i):
			if G[i][j]:
				total_edges += 1

	for b in range(1, 2**N - 1):
		count = 0
		size = 0

		for i in range(N):
			if (powers[i] & b):
				size += 1

				for j in range(i):
					if (powers[j] & b) and G[i][j]:
						count += 1

		metric = max(metric, (count - (size * total_edges) / N) / (size / N * (1.0 - size / N)))
		# metric = max(metric, (count - (size * total_edges) / N))

	print(metric)
	return metric

def list_of_increases(G):
	global N, powers
	G_metric = calculate_metric(G)
	result = []

	for i in range(N):
		best_metric = 100000000.0

		for j in range(N):
			if G[i][j]:
				G[i][j] = False
				G[j][i] = False
				best_metric = min(best_metric, calculate_metric(G))
				G[i][j] = True
				G[j][i] = True

		result.append(best_metric - G_metric)

	return result

while True:
	G = generate_random()

	G_WORKS = True
	for i in range(N):
		WORKS = False
		for j in range(N):
			if G[i][j]:
				WORKS = True
		if not WORKS:
			G_WORKS = False
	if G_WORKS:
		break

print(numpy.array(G))
l = list_of_increases(G)
sum = 0.0
for i in range(N):
	sum += l[i]
print(l)
print(sum)

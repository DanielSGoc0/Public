N = 16
K = 2

E = [[False for j in range(N)] for i in range(N)]
for i in range(N):
	for j in range(N):
		if (N - i + j) % N < K:
			E[i][j] = True

powers = [2**i for i in range(N)]

p = [0 for b in range(2**N)]
p[0] = 1

for b in range(1, 2**N):
	if b % 1024 == 0:
		print(b)
	sum = 0

	for i in range(N):
		best = 0

		for j in range(N):
			if E[i][j] and (b & powers[j]):
				best = max(best, p[b - powers[j]])
		sum += best

	p[b] = sum

print(p[2**N - 1])
print(p[2**N - 1] / N**N)



# N = 2		K = 1
# 2
# 0.5

# N = 4		K = 2
# 112
# 0.4375

# N = 6		K = 3
# 19512
# 0.4182098765432099

# N = 8		K = 4
# 6967296
# 0.415283203125

# N = 10	K = 5
# 4135210150
# 0.413521015

# N = 12	K = 6
# 3677802801984
# 0.4124900592280096

# N = 14	K = 7
# 4590652874294256
# 0.4131254548670353

# N = 16	K = 8
# 7609384593732857344
# 0.4125055653901446

# N = 18	K = 9
# 16265222623744568949744
# 0.4133851962450574

# N = 20	K = 10
# 43310105411817601168996800
# 0.4130373517209778



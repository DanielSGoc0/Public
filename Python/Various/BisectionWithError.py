
N = 160

dp = [[[100000 for c in range(N)] for b in range(N)] for a in range(N)]

dp[1][0][0] = 0
dp[0][1][0] = 0
dp[0][0][1] = 0

for s in range(2, N):
	for b in range(s + 1):
		for a in range(s - b + 1):
			c = s - a - b

			d = 1000000
			for i in range(1, a + 1):
				d = min(d, 1 + max(dp[i + b][0][0], dp[a - i][b][c]))
			for i in range(1, b):
				d = min(d, 1 + max(dp[a][i][b - i], dp[i][b - i][c]))
			for i in range(1, c + 1):
				d = min(d, 1 + max(dp[0][0][i + b], dp[a][b][c - i]))
			
			dp[a][b][c] = d

			# print(a, b, c, d)

for i in range(N):
	print(i, '\t', dp[0][i][0])

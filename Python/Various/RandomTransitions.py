import random

N = 10

LEFT = [random.randint(0, N - 1) for i in range(N)] 
RIGHT = [random.randint(0, N - 1) for i in range(N)] 

for k in range(100):
	NEXT_LEFT = [0 for i in range(N)]
	NEXT_RIGHT = [0 for i in range(N)]

	for i in range(N):
		if random.randint(0, 1) == 0:
			NEXT_LEFT[i] = LEFT[LEFT[i]]
		else:
			# NEXT_LEFT[i] = RIGHT[LEFT[i]]
			NEXT_LEFT[i] = LEFT[RIGHT[i]]
		if random.randint(0, 1) == 0:
			# NEXT_RIGHT[i] = LEFT[RIGHT[i]]
			NEXT_RIGHT[i] = RIGHT[LEFT[i]]
		else:
			NEXT_RIGHT[i] = RIGHT[RIGHT[i]]

	for i in range(N):
		LEFT[i] = NEXT_LEFT[i]
		RIGHT[i] = NEXT_RIGHT[i]

	print(k)
	print(LEFT)
	print(RIGHT)

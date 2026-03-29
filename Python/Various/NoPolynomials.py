# Check if there exists a sequence a_i such that no polynomial
# Q(t) of degree d intersects a in more than d + 1 points.
import pprint

P = 7

inverses = [1 for j in range(P)]
for j in range(1, P):
	for k in range(P - 2):
		inverses[j] = (j * inverses[j]) % P

intersections = [[0 for i in range(P)] for j in range(P)]
a = [i for i in range(P)]

# WLOG a[0] = 0, a[1] = 1
for j in range(2, P):
	intersections[0][j] += 1
	intersections[1][j] += 1
	intersections[j][j] += 1

def get_inv(x):
	global inverses

	if x < 0:
		return -inverses[(-x) % P]
	return inverses[x % P]

def change_intersections(j, diff):
	global P, a, intersections

	# iterate through all non-empty bitmasks of {0, ... , j - 1}
	bitmask = [False for x in range(j + 1)]
	bitmask[j] = True

	while(True):
		# Extrapolate polynomial on points from bitmask and j.
		for j2 in range(j + 1, P):
			# Calculate the extrapolation polynomial at j2.
			sum = 0

			for x in range(j + 1):
				if not bitmask[x]:
					continue

				y = a[x]
				for x2 in range(j + 1):
					if bitmask[x2] and x != x2:
						y = (y * (j2 - x2) * get_inv(x - x2)) % P
				# print(sum, y)
				sum += y

			# Update the intersections array
			sum %= P
			intersections[sum][j2] += diff


		# Generate next bitmask or abort
		END = True
		for x in range(j):
			if not bitmask[x]:
				bitmask[x] = True
				END = False
				break
			bitmask[x] = False
		if END:
			break


def recur(i, j):
	global P, a, intersections

	if j == P:
		print("SUCCESS!")
		for k in range(P):
			print(a[k])
		pprint.pprint(intersections)
		return
	
	if intersections[i][j] > 0:
		return
	
	# Now we assume that we can set a[j] = i
	a[j] = i

	# We remove all subsequent intersections.
	change_intersections(j, 1)

	# We recurse
	for i2 in range(P):
		recur(i2, j + 1)

	# Finally, we remove the intersections from before
	change_intersections(j, -1)


for i in range(P):
	recur(i, 2)

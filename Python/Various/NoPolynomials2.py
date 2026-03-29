# Check if there exists a sequence a_i such that no polynomial
# Q(t) of degree d intersects a in more than d + 1 points.
# However, index 0 is left out (i.e. a_0 does not exist)
# Then a_i = 1/i actually works! But does anything else?
import pprint

P = 19

inverses = [1 for j in range(P)]
for j in range(1, P):
	for k in range(P - 2):
		inverses[j] = (j * inverses[j]) % P

intersections = [[0 for i in range(P)] for j in range(P)]
a = [i for i in range(P)]

# WLOG a[1] = 1, a[2] = 2
a[1] = 1
a[2] = inverses[2]
for j in range(3, P):
	intersections[1][j] += 1
	intersections[inverses[2]][j] += 1
	intersections[(1 + (j - 1) * (inverses[2] - 1)) % P][j] += 1

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
		for x in range(1, j):
			if not bitmask[x]:
				bitmask[x] = True
				END = False
				break
			bitmask[x] = False
		if END:
			break


def recur(i, j):
	global P, a, intersections
	
	if intersections[i][j] > 0:
		return
	
	# Now we assume that we can set a[j] = i
	a[j] = i

	# We remove all subsequent intersections.
	change_intersections(j, 1)

	# We recurse
	if j + 1 == P:
		print("SUCCESS!")
		for k in range(P):
			print(a[k])
		pprint.pprint(intersections)
		return
	else:
		for i2 in range(P):
			recur(i2, j + 1)

	# Finally, we remove the intersections from before
	change_intersections(j, -1)


for i in range(P):
	recur(i, 3)

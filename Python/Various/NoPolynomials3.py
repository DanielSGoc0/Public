# Check if there exists a sequence a_i such that no polynomial
# Q(t) of degree d intersects a in more than d + 1 points.
# However, indices 0 and 1 are left out (i.e. a_0 and a_1 do not exist)
import pprint

P = 17

inverses = [1 for j in range(P)]
for j in range(1, P):
	for k in range(P - 2):
		inverses[j] = (j * inverses[j]) % P

intersections = [[0 for i in range(P)] for j in range(P)]
a = [i for i in range(P)]

# WLOG a[1] = 1, a[2] = 2
a[2] = inverses[2]
a[3] = inverses[3]
for j in range(4, P):
	intersections[inverses[2]][j] += 1
	intersections[inverses[3]][j] += 1
	intersections[(P + inverses[2] + (j - 2) * (inverses[3] - inverses[2])) % P][j] += 1

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
		for x in range(2, j):
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

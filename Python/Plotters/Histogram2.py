# A histogram

import matplotlib.pyplot as plt

N = 512
INTERVAL = 12
LEVELS = 12
filename = "samples/N" + str(N) + "INTERVAL" + str(INTERVAL) + "LEVELS" + str(LEVELS) + ".txt"
input_file = open(filename, 'r')

# ======================================================================

Y = []
DY = []
DDY = []
COUNT = 0

for line in input_file.readlines():
	COUNT += 1
	[y, dy, ddy] = [float(v) for v in line.split(' ')]
	Y.append(y)
	DY.append(dy)
	DDY.append(ddy)

mean = 0.0
deviation = 0.0
for y in Y:
	mean += y
	deviation += y*y
meandy = 0.0
deviationdy = 0.0
for dy in DY:
	meandy += dy
	deviationdy += dy*dy

print(COUNT)
print(mean/COUNT)
print(deviation/COUNT - (mean/COUNT)**2)
print(meandy/COUNT)
print(deviationdy/COUNT - (meandy/COUNT)**2)

# ======================================================================

n, bins, patches = plt.hist(x=Y, bins=100, color='#0504aa', alpha=0.7, rwidth=0.85)
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

input_file.close()

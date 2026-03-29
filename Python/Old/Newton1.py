# Some stuff with Newton differences. Unfinished.

from tkinter import *
import tkinter
import copy
import random
import math


RozmiarX = 1800
RozmiarY = 1000

MoveX = RozmiarX/2
MoveY = RozmiarY/2
Zoom = 0

root = Tk()
frame = Frame(root)
frame.pack()

var1 = StringVar()
var2 = StringVar()
var3 = StringVar()

C = tkinter.Canvas(root, height=RozmiarY, width=RozmiarX, bg = 'white')
C.pack()
TEST = False

NODES = []
N = 0
H = 1


# ===========================================================================================================================
# ===========================================================================================================================
# ===========================================================================================================================


# nodes:
# a list of pairs (x, ((y0, true/false), (y1, true/false), ... , (y(k-1), true/false)))
# each true/false says whether the cooridnate is being kept constant or not. (true if fixed)

# offset: always at least 1

def Newton1(nodes):
	N = 0
	x = []

	for (xi, yi) in nodes:
		N = N + len(yi)
		for k in range(len(yi)):
			x.append(xi)

	A = [[0 for _ in range(N)] for _ in range(N)]
	done = [[False for _ in range(N)] for _ in range(N)]
	
	s = 0
	for (xi, yi) in nodes:
		k = len(yi)
		for a in range(s, s + k):
			for b in range(a, s + k):
				A[a][b] = yi[b - a][0]
				done[a][b] = True
		s = s + k

	for l in range(1, N):
		for a in range(0, N - l):
			b = a + l
			if not done[a][b]:
				A[a][b] = (A[a][b - 1] - A[a + 1][b])/(x[a] - x[b])

	return A


def Newton2(nodes):
	N = 0
	x = []

	for (xi, yi) in nodes:
		N = N + len(yi)
		for k in range(len(yi)):
			x.append(xi)

	A = [[0 for _ in range(N)] for _ in range(N)]
	dA = [[[0 for _ in range(N)] for _ in range(N)] for _ in range(N)]
	done = [[False for _ in range(N)] for _ in range(N)]
	
	s = 0
	for (xi, yi) in nodes:
		k = len(yi)
		for a in range(s, s + k):
			for b in range(a, s + k):
				A[a][b] = yi[b - a][0]
				dA[a][b][s + b - a] = 1
				done[a][b] = True
		s = s + k

	for l in range(1, N):
		for a in range(0, N - l):
			b = a + l
			if not done[a][b]:
				A[a][b] = (A[a][b - 1] - A[a + 1][b])/(x[a] - x[b])
				for i in range(0, N):
					dA[a][b][i] = (dA[a][b - 1][i] - dA[a + 1][b][i])/(x[a] - x[b])

	return (A, dA)


# I assume that offset >= 1.
def norm_sup1(nodes, offset):
	A = Newton1(nodes)
	N = len(A)

	sup = 0
	for k in range(offset, N):
		p = 1.0/k

		for a in range(0, N - k):
			sup = max(sup, abs(A[a][a + k]) ** p)

	return sup


# I assume that offset >= 1.
def norm_sup2(nodes, offset):
	(A, dA) = Newton2(nodes)
	N = len(A)

	sup = 0
	dsup = [0 for _ in range(N)]
	for k in range(offset, N):
		p = 1.0/k

		for a in range(0, N - k):
			val = abs(A[a][a + k]) ** p
			sign = 1

			if A[a][a + k] < 0:
				sign = -1

			if val > sup:
				sup = val

				for i in range(N):
					dsup[i] = p * dA[a][a + k][i] * val / abs(A[a][a + k]) * sign

	return (sup, dsup)


def norm_l1_sup1(nodes, offset):
	A = Newton1(nodes)
	N = len(A)

	sum = 0.0

	for k in range(offset, N):
		p = 1.0/k
		sup = 0.0

		for a in range(0, N - k):
			sup = max(sup, abs(A[a][a + k]) ** p)

		sum = sum + sup

	return sum


def norm_l1_sup2(nodes, offset):
	(A, dA) = Newton2(nodes)
	N = len(A)

	sum = 0.0
	dsum = [0.0 for _ in range(N)]

	for k in range(offset, N):
		p = 1.0/k
		sup = 0.0
		dsup = [0.0 for _ in range(N)]

		for a in range(0, N - k):
			val = abs(A[a][a + k]) ** p
			sign = 1

			if A[a][a + k] < 0:
				sign = -1

			if val > sup:
				sup = val

				for i in range(N):
					dsup[i] = p * dA[a][a + k][i] * val / abs(A[a][a + k]) * sign

		sum = sum + sup
		for i in range(N):
			dsum[i] = dsum[i] + dsup[i]

	return (sum, dsum)


def upgrade(nodes, offset, h):
	(norm, dnorm) = norm_l1_sup2(nodes, offset)
	N = len(dnorm)

	# a simple gradient descent. Along every coordinate.
	SUM = 0.0

	i = 0
	j = 0
	k = 0
	while k < N:
		if not nodes[i][1][j][1]:
			SUM = SUM + dnorm[i]**2

		j = j + 1
		if j == len(nodes[i][1]):
			j = 0
			i = i + 1
		k = k + 1

	SUM = SUM ** (1.0/2)

	if SUM < 0.001:
		return (nodes, norm, False)
	
	
	upgraded = copy.deepcopy(nodes)

	i = 0
	j = 0
	k = 0
	while k < N:
		if not nodes[i][1][j][1]:
			upgraded[i][1][j][0] = upgraded[i][1][j][0] - dnorm[i] * h / SUM

		j = j + 1
		if j == len(nodes[i][1]):
			j = 0
			i = i + 1
		k = k + 1

	return (upgraded, norm, True)



def F(x):
	return 1/(1+x**2)

def create_random(n):
	x = []

	for _ in range(n):
		x.append(random.gauss(0, 1))
	x.sort()

	nodes = []

	for xi in x:
		yi = F(xi)
		b = random.getrandbits(1)
		if b == 1:
			nodes.append([xi, [[yi, b]]])
		else:
			nodes.append([xi, [[yi + random.gauss(0, 0.1), b]]])

	return nodes


# ===========================================================================================================================
# ===========================================================================================================================
# ===========================================================================================================================


def test():
	global MoveX, MoveY, X_axis, Y_axis, RozmiarX, RozmiarY, NODES, N
	nodes = []
	K = 100

	for i in range(K):
		x = (1.0 + i)/(1.0 + K)
		y = math.sin(1.0/x)
		# y = 1.0/(1.0 + (x-1.0/2)*(x-1.0/2))
		# y = math.exp(x)
		# if i == K/2:
			# y = y + 0.00001
		nodes.append([x, [(y, True)]])

	results = Newton1(nodes)
	# print(results)

	C.delete(ALL)

	X_axis = C.create_line(0, MoveY, RozmiarX, MoveY, fill = 'red')
	Y_axis = C.create_line(MoveX, 0, MoveX, RozmiarY, fill = 'red')

	for i in range(K):
		for j in range(i, K):
			if j - i == N:
				#PX1 = nodes[i][0] * 2**(Zoom) + MoveX
				#PX2 = nodes[j][0] * 2**(Zoom) + MoveX
				PX = (nodes[i][0] + nodes[j][0])/2.0 * 2**(Zoom) + MoveX
				
				PY = results[i][j]
				if PY >= 0:
					PY = abs(PY)**(1.0/N)
				else:
					PY = -abs(PY)**(1.0/N)
				PY = MoveY - PY * 2**Zoom
				#C.create_line(PX1 - 1, PY, PX2 + 1, PY, fill = 'red', width=5)
				C.create_oval(PX - 1, PY, PX + 1, PY, fill = 'red', width=5)




def Action1():
	global MoveX, MoveY
	MoveX -= 100
	paint()

def Action2():
	global MoveX, MoveY
	MoveX += 100
	paint()

def Action3():
	global MoveX, MoveY
	MoveY += 100
	paint()

def Action4():
	global MoveX, MoveY
	MoveY -= 100
	paint()

def Action5():
	global Zoom, var2
	Zoom += 1
	var2.set("Zoom = " + str(Zoom))
	paint()

def Action6():
	global Zoom, var2
	Zoom -= 1
	var2.set("Zoom = " + str(Zoom))
	paint()

def Action7():
	global H, var3
	H = H * 2
	var3.set("H = " + str(H))
	paint()

def Action8():
	global H, var3
	H = H / 2
	var3.set("H = " + str(H))
	paint()

def Action9():
	global N, var1
	N += 1
	var1.set("N = " + str(N))
	paint()

def Action10():
	global N, var1
	if N > 0:
		N -= 1
		var1.set("N = " + str(N))
	paint()

def Action11():
	global NODES
	NODES = create_random(N)
	paint()

def Action12():
	global NODES, H
	NODES = upgrade(NODES, 1, H)[0]
	paint()

def Action13():
	global TEST
	if TEST:
		TEST = False
	else:
		TEST = True
	paint()


def paint():
	global MoveX, MoveY, X_axis, Y_axis, RozmiarX, RozmiarY, NODES
	if TEST:
		test()
	else:
		C.delete(ALL)

		X_axis = C.create_line(0, MoveY, RozmiarX, MoveY, fill = 'red')
		Y_axis = C.create_line(MoveX, 0, MoveX, RozmiarY, fill = 'red')

		for (xi, yi) in NODES:
			PX = xi * 2**(Zoom) + MoveX
			PY = MoveY - yi[0][0] * 2**Zoom
			if yi[0][1]:
				C.create_oval(PX - 2, PY - 2, PX + 2, PY + 2, fill = 'black')
			else:
				C.create_oval(PX - 2, PY - 2, PX + 2, PY + 2, fill = 'white')



Button1 = tkinter.Button(frame, text ="RIGHT", command = Action1)
Button1.pack(side = LEFT)

Button2 = tkinter.Button(frame, text ="LEFT", command = Action2)
Button2.pack(side = LEFT)

Button3 = tkinter.Button(frame, text ="UP", command = Action3)
Button3.pack(side = LEFT)

Button4 = tkinter.Button(frame, text ="DOWN", command = Action4)
Button4.pack(side = LEFT)

Button5 = tkinter.Button(frame, text ="IN", command = Action5)
Button5.pack(side = LEFT)

Button6 = tkinter.Button(frame, text ="OUT", command = Action6)
Button6.pack(side = LEFT)

Button7 = tkinter.Button(frame, text ="H+", command = Action7)
Button7.pack(side = LEFT)

Button8 = tkinter.Button(frame, text ="H-", command = Action8)
Button8.pack(side = LEFT)

Button9 = tkinter.Button(frame, text ="N++", command = Action9)
Button9.pack(side = LEFT)

Button10 = tkinter.Button(frame, text ="N--", command = Action10)
Button10.pack(side = LEFT)

Button11 = tkinter.Button(frame, text ="RANDOM", command = Action11)
Button11.pack(side = LEFT)

Button12 = tkinter.Button(frame, text ="UPGRADE", command = Action12)
Button12.pack(side = LEFT)

Button13 = tkinter.Button(frame, text ="TEST", command = Action13)
Button13.pack(side = LEFT)


Label1 = Label(frame, textvariable=var1, relief=RAISED)
var1.set("N = 0")
Label1.pack(side = LEFT)

Label2 = Label(frame, textvariable=var2, relief=RAISED)
var2.set("Zoom = 0")
Label2.pack(side = LEFT)

Label3 = Label(frame, textvariable=var3, relief=RAISED)
var3.set("H = " + str(H))
Label3.pack(side = LEFT)

root.mainloop()




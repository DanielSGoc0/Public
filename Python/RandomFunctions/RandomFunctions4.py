from tkinter import *
import tkinter
import copy
import random
import math
import numpy
import scipy

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

N = 16
coeffs = [0 for i in range(N)]
SOLS = []

# lnfactorials = [scipy.special.gammaln(n + 1) for n in range(1024)]
sqrts = [(n + 1) ** (-1/2.0) for n in range(16384)]


def new_random_function():
	global N, coeffs

	coeffs = numpy.random.normal(0, 1, N)
	for i in range(N):
		if i % 2 == 1:
			coeffs[i] = 0
		else:
			coeffs[i] = abs(coeffs[i])

# def f(x):
# 	global N, coeffs
# 	y = 0
# 	for i in range(N):
# 		z = 0
# 		if x == 0:
# 			if i == 0:
# 				z = 1
# 			else:
# 				z = 0
# 		else:
# 			z = - x*x / 2.0
# 			z = z - lnfactorials[i] / 2.0
# 			z = z + i * math.log(abs(x))
# 			z = math.exp(z)
# 			if x < 0 and i % 2 == 1:
# 				z = -z
# 		y = y + z * coeffs[i]
# 	return y

def f(x):
	global N, coeffs
	y = 0
	for i in range(N - 1, -1, -1):
		y = y * x * sqrts[i]
		y = y + coeffs[i]
	return y * math.exp(- x*x/2.0)
	

def map_point(x, y):
	global MoveX, MoveY, Zoom
	x = MoveX + x * 2**Zoom
	y = MoveY - y * 2**Zoom
	return (x, y)

def evaluate(x):
	global MoveX, MoveY, Zoom

	if x % 2 == 0:
		return x
	else:
		x -= MoveX
		x *= 2**(-Zoom)

		y = f(x)

		return MoveY - y * 2**Zoom

def flow_equation(T):
	x = 0
	for _ in range(1000):
		dx = f(x) * T/1000
		x = x + dx
	return x


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
	global SOLS
	# for i in range(100):
	new_random_function()

	# X = flow_equation(50)
	# Y = f(X)
	# SOLS.append((X, Y))

	paint()
	
def Action9():
	global N, var1
	N = N * 2
	var1.set("N = " + str(N))
	# paint()

def Action10():
	global N, var1
	if N > 1:
		N = N//2
		var1.set("N = " + str(N))
	# paint()



def paint():
	global MoveX, MoveY, line_new, cord, X_axis, Y_axis, RozmiarX, RozmiarY, SOLS
	C.delete(ALL)

	for w in range(0, 160):
		cord = [evaluate(i) for i in range(w*12, (w + 1)*12 + 2)]
		line_new = C.create_line(cord)
	#print(MoveY - funkcja(11.5 + 0.3333333 + MoveX))
	X_axis = C.create_line(0, MoveY, RozmiarX, MoveY, fill = 'red')
	Y_axis = C.create_line(MoveX, 0, MoveX, RozmiarY, fill = 'red')

	# for (X, Y) in SOLS:
	# 	(X, Y) = map_point(X, Y)
	# 	C.create_oval(X - 1, Y - 1, X + 1, Y + 1, fill = 'red')
	# print(SOLS)

# def paint():
# 	global MoveX, MoveY, X_axis, Y_axis, RozmiarX, RozmiarY, NODES
# 	C.delete(ALL)

# 	X_axis = C.create_line(0, MoveY, RozmiarX, MoveY, fill = 'red')
# 	Y_axis = C.create_line(MoveX, 0, MoveX, RozmiarY, fill = 'red')

# 	for i in range(len(NODES) - 1):
# 		PX0 = NODES[i][0] * 2**(Zoom) + MoveX
# 		PY0 = MoveY - NODES[i][1] * 2**Zoom
# 		PX1 = NODES[i + 1][0] * 2**(Zoom) + MoveX
# 		PY1 = MoveY - NODES[i + 1][1] * 2**Zoom

# 		C.create_line(PX0, PY0, PX1, PY1, fill = 'black')



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

Button7 = tkinter.Button(frame, text ="RANDOM", command = Action7)
Button7.pack(side = LEFT)

Button9 = tkinter.Button(frame, text ="N++", command = Action9)
Button9.pack(side = LEFT)

Button10 = tkinter.Button(frame, text ="N--", command = Action10)
Button10.pack(side = LEFT)



Label1 = Label(frame, textvariable=var1, relief=RAISED)
var1.set("N = 16")
Label1.pack(side = LEFT)

Label2 = Label(frame, textvariable=var2, relief=RAISED)
var2.set("Zoom = 0")
Label2.pack(side = LEFT)

root.mainloop()




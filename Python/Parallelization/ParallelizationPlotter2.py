# Plots function F, as target point changes

from tkinter import *
import tkinter
import math
import numpy as np
from numpy import *
import scipy.linalg


RozmiarX = 1800
RozmiarY = 1000

MoveX = 0.0
MoveY = 0.0
Zoom = 8

root = Tk()
frame = Frame(root)
frame.pack()

var1 = StringVar()

C = tkinter.Canvas(root, height=RozmiarY, width=RozmiarX, bg = 'white')
C.pack()
N = 64
NODES = []
POINTS = []

def K(x, y):
	return np.exp(-np.dot(x - y, x - y) / 2.0)

def f(X, Z, err2, X_target):
	# Z is a matrix of dimensions D x N.
	S = len(X)
	D = len(X_target)

	Sigma = np.zeros((S+1, S+1), dtype='double')

	for i in range(S):
		for j in range(S):
			Sigma[i][j] = K(X[i], X[j])
		Sigma[i][i] += err2[i]

		Sigma[i][S] = K(X[i], X_target)
		Sigma[S][i] = Sigma[i][S]
	Sigma[S][S] = K(X_target, X_target) + 1.0

	Chol = scipy.linalg.cholesky(Sigma, lower=True)
	
	ANS = 0.0
	for i in range(S):
		for j in range(D):
			ANS += (X_target[j] - X[i][j]) * Z[j][i] * Chol[S][i]

	return ANS

def new_plot():
	global N, NODES, POINTS
	NODES.clear()
	POINTS.clear()

	X = np.array([[-0.3], [0.3]])
	Z = np.array([[scale1.get()/1000.0, scale2.get()/1000.0]])
	err2 = [0.0, 0.0]

	for k in range(-10*N, 10*N + 1):
		x = np.array([(k * 1.0)/N])
		val = f(X, Z, err2, x)
		NODES.append((x[0], val))

	POINTS.append((0.3, f(X, Z, err2, np.array([0.3]))))
	POINTS.append((-0.3, f(X, Z, err2, np.array([-0.3]))))


def Action1():
	global MoveX, MoveY, Zoom
	MoveX += 100 * 2**(-Zoom)
	paint()

def Action2():
	global MoveX, MoveY, Zoom
	MoveX -= 100 * 2**(-Zoom)
	paint()

def Action3():
	global MoveX, MoveY, Zoom
	MoveY += 100 * 2**(-Zoom)
	paint()

def Action4():
	global MoveX, MoveY, Zoom
	MoveY -= 100 * 2**(-Zoom)
	paint()

def Action5():
	global Zoom, var1
	Zoom += 1
	var1.set("Zoom = " + str(Zoom))
	paint()

def Action6():
	global Zoom, var1
	Zoom -= 1
	var1.set("Zoom = " + str(Zoom))
	paint()

def Action7():
	new_plot()
	paint()


def slider(t):
	new_plot()
	paint()

def paint_coords(x, y):
	global MoveX, MoveY, RozmiarX, RozmiarY, Zoom
	return ((x - MoveX) * 2**Zoom + RozmiarX/2, RozmiarY/2 - (y - MoveY) * 2**Zoom)

def paint():
	global MoveX, MoveY, X_axis, Y_axis, RozmiarX, RozmiarY, NODES, POINTS
	C.delete(ALL)

	X_axis = C.create_line(0, MoveY*2**Zoom + RozmiarY/2, RozmiarX, MoveY*2**Zoom + RozmiarY/2, fill = 'red')
	Y_axis = C.create_line(-MoveX*2**Zoom + RozmiarX/2, 0, -MoveX*2**Zoom + RozmiarX/2, RozmiarY, fill = 'red')

	for i in range(len(NODES) - 1):
		(PX0, PY0) = paint_coords(NODES[i][0], NODES[i][1])
		(PX1, PY1) = paint_coords(NODES[i + 1][0], NODES[i + 1][1])

		C.create_line(PX0, PY0, PX1, PY1, fill = 'black')

	for (x, y) in POINTS:
		(px, py) = paint_coords(x, y)
		C.create_oval(px - 2, py - 2, px + 2, py + 2)


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

Button7 = tkinter.Button(frame, text ="NEW PLOT", command = Action7)
Button7.pack(side = LEFT)


Label1 = Label(frame, textvariable=var1, relief=RAISED)
var1.set("Zoom = " + str(Zoom))
Label1.pack(side = LEFT)

scale1 = Scale(frame, from_=-200, to=200, orient=HORIZONTAL, command=slider)
scale1.set(0)
scale1.pack(side = LEFT)

scale2 = Scale(frame, from_=-200, to=200, orient=HORIZONTAL, command=slider)
scale2.set(0)
scale2.pack(side = LEFT)

scale3 = Scale(frame, from_=-200, to=200, orient=HORIZONTAL, command=slider)
scale3.set(0)
scale3.pack(side = LEFT)

root.mainloop()

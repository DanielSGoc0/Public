# This version of plotter has variable sliders and better zoom!
# So basically we just condition Gaussian and plot condition covariance wrt fixed point.

from tkinter import *
import tkinter
import math
import numpy as np
from numpy import *


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
NODES = []
POINTS = []

# ===========================================================================================

# we take conditioned covariance of points n and n+1
def calc(X, p):
	n = len(p)

	Sigma = np.zeros((n + 2, n + 2), dtype='double')
	for i in range(n + 2):
		for j in range(n + 2):
			Sigma[i][j] = np.exp(-np.dot(X[i] - X[j], X[i] - X[j]) / 2.0)
	for i in range(n):
		Sigma[i][i] += 1.0 / p[i]
	C = np.linalg.cholesky(Sigma)

	return C[n + 1][n] * C[n][n]


def new_plot():
	global NODES, POINTS, scale_x1, scale_p1, scale_x2, scale_p2
	NODES.clear()
	POINTS.clear()

	x1 = scale_x1.get() / 100.0
	p1 = scale_p1.get() / 100.0
	x2 = scale_x2.get() / 100.0
	p2 = scale_p2.get() / 100.0

	LINE = []
	LINE2 = []
	LINE3 = []
	for scale_x3 in range(-300, 301):
		x = scale_x3 / 100.0
		x += 0.005
		X = np.array([[x1], [x2], [x], [0.0]], dtype='double')
		p = np.array([p1, p2], dtype='double')
		LINE.append((x, calc(X, p)))
		LINE2.append((x, np.exp(-np.dot(X[2] - X[3], X[2] - X[3]) / 2.0)))
		LINE3.append((x, -np.exp(-np.dot(X[2] - X[3], X[2] - X[3]) / 2.0)))

	NODES.append(LINE)
	NODES.append(LINE2)
	NODES.append(LINE3)


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

	for LINE in NODES:
		for i in range(len(LINE) - 1):
			(PX0, PY0) = paint_coords(LINE[i][0], LINE[i][1])
			(PX1, PY1) = paint_coords(LINE[i + 1][0], LINE[i + 1][1])

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

scale_x1 = Scale(frame, from_=-200, to=200, orient=HORIZONTAL, command=slider)
scale_x1.set(0)
scale_x1.pack(side = LEFT)

scale_p1 = Scale(frame, from_=1, to=500, orient=HORIZONTAL, command=slider)
scale_p1.set(0)
scale_p1.pack(side = LEFT)

scale_x2 = Scale(frame, from_=-200, to=200, orient=HORIZONTAL, command=slider)
scale_x2.set(0)
scale_x2.pack(side = LEFT)

scale_p2 = Scale(frame, from_=1, to=500, orient=HORIZONTAL, command=slider)
scale_p2.set(0)
scale_p2.pack(side = LEFT)


def key_handler(event):
	if event.keycode == 114:
		Action1()
	elif event.keycode == 113:
		Action2()
	elif event.keycode == 111:
		Action3()
	elif event.keycode == 116:
		Action4()
	elif event.keycode == 31:
		Action5()
	elif event.keycode == 32:
		Action6()

root.bind("<Key>", key_handler)

root.mainloop()

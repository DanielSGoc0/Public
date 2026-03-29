# This version of plotter has variable sliders and better zoom!

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

def phi(r2):
	# return np.exp(-r2) + 0.1 * np.exp(-r2 / 4.0)
	return np.exp(-r2 / 2.0)

def calc1(X, sigma2, l):
	s = len(l) - 1
	n = l[0]
	N = l[s]

	Sigma = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N + 1):
		for j in range(N + 1):
			Sigma[i][j] = phi(np.dot(X[i] - X[j], X[i] - X[j]))
		Sigma[i][i] += sigma2[i]
	C = np.linalg.cholesky(Sigma)

	# start by calculating ANS
	ANS = 0.0
	for k in range(s):
		ANS += X[N][n + k] * np.sqrt(np.dot(C[N, l[k]:l[k+1]], C[N, l[k]:l[k+1]]))
	
	# Next calculate weigths w
	w = np.zeros((s), dtype='double')
	for k in range(s - 1):
		w[k] += X[N][n + k] / np.sqrt(np.dot(C[N, l[k]:l[k+1]], C[N, l[k]:l[k+1]]))
		w[k] -= X[N][n + k + 1] / np.sqrt(np.dot(C[N, l[k+1]:l[k+2]], C[N, l[k+1]:l[k+2]]))
	w[s-1] = X[N][n + s-1] / np.sqrt(np.dot(C[N, l[s-1]:l[s]], C[N, l[s-1]:l[s]]))

	return (ANS, w)


def calc2(X, sigma2, l, ANS, w):
	s = len(l) - 1
	n = l[0]
	N = l[s]

	Sigma = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N + 1):
		for j in range(N + 1):
			Sigma[i][j] = phi(np.dot(X[i] - X[j], X[i] - X[j]))
		Sigma[i][i] += sigma2[i]
	C = np.linalg.cholesky(Sigma)

	ANS2 = 0.0
	for k in range(s):
		ANS2 += w[k] * np.dot(C[N, n:l[k+1]], C[N, n:l[k+1]])
	ANS2 = (ANS + ANS2) / 2.0

	return ANS2


def get_X(x):
	# X0 = np.array([[1.0, 0.0, 0.0, 0.0], [-1.0, 0.1, 0.0, 0.0], [0.8, -0.5, 0.0, 0.0], [x, -0.3, 0.5+x, 0.0], [0.0, 0.0, 0.3, 1.0]])
	X0 = np.array([[1.0, 0.0, 0.0, 0.0], [-1.0, 0.1, 0.0, 0.0], [0.3, -0.5 - x*0.05, 0.0, 0.0], [x, -0.3, 0.5+x, 0.0], [0.0, 0.0, 0.3, 0.6]])
	# X0 = np.array([[1.0, 0.0, 0.0, 0.0], [-1.0, 0.1, 0.0, 0.0], [0.8, -0.5, 0.0, 0.0], [x, -0.3, 0.5+x, 0.0], [0.0, 0.0, 0.3 + 0.3*x, np.sqrt(1.0 - (0.3 + 0.3*x)**2)]])
	# X0 = np.array([[1.0, 0.0, 0.0, 0.0], [-1.0, 0.1, 0.0, 0.0], [0.8, -0.5, 0.0, 0.0], [0.0, -0.3, 0.5, 0.0], [0.0, 0.0, 0.3 + 0.3*x, np.sqrt(1.0 - (0.3 + 0.3*x)**2)]])
	return X0

def get_sigma2(x):
	# sigma2 = np.array([0.0, 0.0, 1.0, 2.0, 0.0])
	# sigma2 = np.array([0.0, 0.0, 1.0 + x**2, 2.0 + 0.3*x**2, 0.0])
	sigma2 = np.array([(x + 0.5)**2 * 0.1, (x - 0.5)**2 * 0.07, 1.0 + x**2, 2.0 + 0.3*x**2, 0.0])
	return sigma2

def new_plot():
	global NODES, POINTS, scale_A, scale_F, scale_R
	NODES.clear()
	POINTS.clear()

	l = [2, 3, 4]
	F = scale_F.get() / 100.0
	X0 = get_X(F)
	sigma0 = get_sigma2(F)

	(ANS, w) = calc1(X0, sigma0, l)
	print(w)

	LINE = []
	for X in range(-100, 101):
		x = X/10.0
		X1 = get_X(x)
		sigma1 = get_sigma2(x)
		(a0, w0) = calc1(X1, sigma1, l)
		w_pos = True
		# for v in w0:
		# 	if v < 0.0:
		# 		w_pos = False
		if w_pos:
			LINE.append((x, 10.0 * a0))
		else:
			LINE.append((x, 0.0))
	NODES.append(LINE)

	LINE = []
	for X in range(-100, 101):
		x = X/10.0
		X1 = get_X(x)
		sigma1 = get_sigma2(x)
		LINE.append((x, 10.0 * calc2(X1, sigma1, l, ANS, w)))
	NODES.append(LINE)



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

scale_F = Scale(frame, from_=-200, to=250, orient=HORIZONTAL, command=slider)
scale_F.set(0)
scale_F.pack(side = LEFT)

scale_A = Scale(frame, from_=0, to=100, orient=HORIZONTAL, command=slider)
scale_A.set(0)
scale_A.pack(side = LEFT)

scale_R = Scale(frame, from_=0, to=100, orient=HORIZONTAL, command=slider)
scale_R.set(0)
scale_R.pack(side = LEFT)


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

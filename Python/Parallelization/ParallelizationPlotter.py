# Plots function G

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
PLOT_G = True
NODES = []
POINTS = []

def K(x, y):
	return np.exp(-np.dot(x - y, x - y)/ 2.0)

def f(X, err2, X_inter):
	# indices 0 to N-1 are data (including initial and evaluation point)
	# index N is the target
	N = len(X) - 1

	Sigma = np.zeros((N+1, N+1), dtype='double')
	Sigma_inter = np.zeros((N+1), dtype='double')
	for i in range(N+1):
		for j in range(N+1):
			Sigma[i][j] = K(X[i], X[j])
		Sigma_inter[i] = K(X[i], X_inter)

	for i in range(N):
		Sigma[i][i] += err2[i]

	return K(X[N], X_inter) - np.dot(Sigma[N][:N], scipy.linalg.solve(Sigma[:N, :N], Sigma_inter[:N], assume_a='pos'))	

def new_plot():
	global N, NODES, POINTS, PLOT_G
	NODES.clear()
	POINTS.clear()

	X = np.array([-0.3, 0.3, scale3.get()/100.0, 0.0])
	err2 = [0.0, 0.0, 1.0]
	
	if PLOT_G:
		maximum_at = 0.0
		maximum = 0.0

		for k in range(-10*N, 10*N + 1):
			x = (k * 1.0)/N

			val = f(X, err2, x)

			if abs(val) > abs(maximum):
				maximum = val
				maximum_at = x

			NODES.append((x, 10*val))

		print(maximum_at)
		POINTS.append((maximum_at, 10.0*maximum))
		POINTS.append((scale3.get()/100.0, 10.0*f(X, err2, scale3.get()/100.0)))
	
	else:
		minimum_at = 0.0
		minimum = 1.0
		for k in range(-10*N, 10*N + 1):
			x = (k * 1.0)/N

			X[2] = x

			val = f(X, err2, 0.0)

			if val < minimum:
				minimum = val
				minimum_at = x

			NODES.append((x, 100.0*val))

		print(minimum_at)
		POINTS.append((minimum_at, 100.0*minimum))



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

def Action8():
	global PLOT_G
	PLOT_G = not PLOT_G
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

Button8 = tkinter.Button(frame, text ="PLOT G", command = Action8)
Button8.pack(side = LEFT)


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

# for k in range(1, 1000):
# 	F_new = k/100.0
# 	for l in range(0, 100):
# 		A_new = (2*l + 1) / 200.0

# 		val0 = next_iter(F_new, A_new)
# 		val1 = next_iter(F_new + 1.0/1000, A_new)
# 		val2 = next_iter(F_new, A_new + 1.0/1000)

# 		w = (val1 - val0)/(val2 - val0) / (2.0 * F_new) - 1.0

# 		# if w < 0.0:
# 		print(F_new, A_new, w)

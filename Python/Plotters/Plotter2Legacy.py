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
N = 64
NODES = []
POINTS = []

F = 0.0
A = 1.0

def f(x, y):
	global F, A

	return (np.exp(-x*x / 2.0) * (F + y * x), 1.0 - np.exp(-x*x) * (1.0 - A + y**2))



def new_plot():
	global N, NODES, POINTS, F, A
	NODES.clear()
	POINTS.clear()

	F = scale_F.get() / 100.0 # F
	A = scale_A.get() / 50.0 # x^2
	
	# S = scale_F.get() / 100.0 * np.exp(1/2.0) # F
	# T = scale_A.get() / 50.0 # x^2

	# F = S * np.exp(-T/2.0) * (1 + T)
	# A = 1.0 - S**2 * np.exp(-T) * T

	# POINTS.append((S * np.exp(-1/2.0) * 2, 1.0 - S**2 * np.exp(-1.0)))

	POINTS.append((F, A))

	# for s0 in range(30):
	# 	S0 = s0/5.0
	# 	LINE = []
	# 	for t0 in range(0, 101):
	# 		T0 = t0/10.0
	# 		LINE.append((S0 * np.exp(-T0/2.0) * (1 + T0), 1.0 - S0**2 * np.exp(-T0) * T0))
	# 	NODES.append(LINE)

	# for t in range(6):
	# 	x = t/10.0
	# 	LINE = []

	# 	for s in range(0, 1001):
	# 		y = s/100.0 * np.sqrt(A)
	# 		# y = s/100.0 * np.sqrt(S**2 * (1 - T) * np.exp(-T))
	# 		try:
	# 			LINE.append(f(x, y))
	# 		except:
	# 			print()
	# 		if y > S:
	# 			break
	# 	NODES.append(LINE)

	# for t in range(-50, 51):
	# 	x = t/10.0
	# 	LINE = []

	# 	for s in range(-10, 11):
	# 		y = s/10.0 * np.sqrt(A)
	# 		# y = s/100.0 * np.sqrt(S**2 * (1 - T) * np.exp(-T))
	# 		try:
	# 			LINE.append(f(x, y))
	# 		except:
	# 			print()
	# 	NODES.append(LINE)

	for C in range(0, 21):
		# c = C/20.0 + 1.0/(C/20.0) - 2.0
		c = tan(C/20.0 * np.pi / 2.0)
		LINE = []

		for s in range(0, 101):
			y = s/100.0 * np.sqrt(A)
			x = c*y
			# y = s/100.0 * np.sqrt(S**2 * (1 - T) * np.exp(-T))
			try:
				LINE.append(f(x, y))
			except:
				print()
		NODES.append(LINE)

	
	# for s in range(0, 3):
	# 	y = s/2.0 * np.sqrt(A)

	# 	LINE = []
	# 	for t in range(-500, 501):
	# 		x = t/100.0
	# 		try:
	# 			LINE.append(f(x, y))
	# 		except:
	# 			print()

	# 	NODES.append(LINE)

	LINE = []
	for s in range(201):
		x = s/100.0
		LINE.append((x, 1.0 - x**2/4))
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

scale_F = Scale(frame, from_=0, to=200, orient=HORIZONTAL, command=slider)
scale_F.set(0)
scale_F.pack(side = LEFT)

scale_A = Scale(frame, from_=0, to=100, orient=HORIZONTAL, command=slider)
scale_A.set(0)
scale_A.pack(side = LEFT)

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

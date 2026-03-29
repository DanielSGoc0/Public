# This version of plotter has variable sliders and better zoom!
# Just another plotter

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

def f(x, y, z):
	x += 0.005
	# x = np.abs(x)
	# return ((x - y) * np.log(x/y) + (y - z) * np.log(y/z)) / ((z - x) * np.log(z/x))
	return 1.0 - (abs(x*y - 1.0/(x*y)) + 1) / ((abs(x - 1.0/x) + 1) * (abs(y - 1.0/y) + 1))

def new_plot():
	global NODES, POINTS, scale_A, scale_F, scale_R
	NODES.clear()
	POINTS.clear()

	y = scale_F.get() / 100.0
	z = scale_A.get() / 100.0

	MIN = 10.0

	LINE = []
	for X in range(-100, 101):
		x = X/50.0
		LINE.append((x, f(x, y, z)))
		MIN = min(LINE[-1][1], MIN)
	print(MIN)
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

scale_F = Scale(frame, from_=-100, to=250, orient=HORIZONTAL, command=slider)
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

#!/usr/bin/python3

from Tkinter import *
import Tkinter
from math import *
from scipy.special import *
from numpy import *
#from bigfloat import *

RozmiarX = 1800
RozmiarY = 1000

MoveX = RozmiarX/2
MoveY = RozmiarY/2

Zoom = 5
N = 0
PI = 3.14159

root = Tk()
frame = Frame(root)
frame.pack()

var1 = StringVar()
var2 = StringVar()
cord = [0, 0, 0, 0]

C = Tkinter.Canvas(root, height=RozmiarY, width=RozmiarX, bg = 'white')
C.pack()

line_old = C.create_line(cord, fill = 'white')
line_new = C.create_line(cord)

X_axis = C.create_line(cord, fill = 'red')
Y_axis = C.create_line(cord, fill = 'red')

DIFF = False
DOTS = True

c = 1
wx = 0.0
wy = 0.0
sign1 = TRUE
sign2 = TRUE
x0 = 0.0
x1 = 0.0
x2 = 0.0
y0 = 0.0
y1 = 0.0
y2 = 0.0

def random_w():
	global wx, wy
	wx = random.uniform(0.5, 2)
	wy = random.uniform(0.5, 2)


def random_initial():
	global x0, y0, sign1, sign2, PI
	a = random.uniform(0, 2*PI)
	b = random.uniform(0, 2*PI)
	x0 = sin(a)*sin(b)
	y0 = cos(a)*cos(b)
	if random.uniform(0, 2) > 1:
		sign1 = FALSE
	if random.uniform(0, 2) > 1:
		sign2 = FALSE

def make_next():
	global x0, y0, x1, y1, sign1, sign2, wx, wy, c
	plus = abs(wx - wy)*sqrt(1 - (x0 + y0)*(x0 + y0))
	minus = abs(wx + wy)*sqrt(1 - (x0 - y0)*(x0 - y0))
	if sign1:
		plus = -plus
	if sign2:
		minus = -minus
	x1 = x0 + c*(plus + minus)/2
	y1 = y0 + c*(plus - minus)/2


def random_next():
	global x1, y1, PI
	a = random.uniform(0, 2*PI)
	b = random.uniform(0, 2*PI)
	x1 = sin(a)*sin(b)
	y1 = cos(a)*cos(b)

def initialize():
	random_w()
	random_initial()
	random_next()

def next_value():
	global wx, wy, c, x0, x1, x2, y0, y1, y2
	S = -(wx*wx + wy*wy)
	T = 2.0*wx*wy
	K = c*c/2.0
	x2 = 2.0*x1 - x0 + K*(S*x1 + T*y1)
	y2 = 2.0*y1 - y0 + K*(S*y1 + T*x1)

	x0 = x1
	x1 = x2
	y0 = y1
	y1 = y2

def funkcja(x):
	global Movex, MoveY, N, seq, x0, y0
	if x % 2 == 0:
		next_value()
		return MoveX - x0 * 2**Zoom
	else:
		return MoveY - y0 * 2**Zoom

def Action1():
	global MoveX, MoveY
	MoveX -= 100

def Action2():
	global MoveX, MoveY
	MoveX += 100

def Action3():
	global MoveX, MoveY
	MoveY += 100

def Action4():
	global MoveX, MoveY
	MoveY -= 100

def Action5():
	global Zoom
	Zoom += 1
	var2.set("Zoom = " + str(Zoom))

def Action6():
	global Zoom
	Zoom -= 1
	var2.set("Zoom = " + str(Zoom))

def Action7():
	global N, var
	if N < 15:
		N += 1
		var1.set("N = " + str(N))

def Action8():
	global N, var
	if N > 0:
		N -= 1
		var1.set("N = " + str(N))

def Action9():
	paint()

def Action10():
	initialize()

def paint():
	global MoveX, MoveY, line_new, cord, X_axis, Y_axis, RozmiarX, RozmiarY, AMT
	C.delete(ALL)
	print(x0)
	print(y0)
	print(x1)
	print(y1)
	for w in range(0, 2**N):
		cord = [funkcja(i) for i in range(0, 18)]
		line_new = C.create_line(cord)
	#print(MoveY - funkcja(11.5 + 0.3333333 + MoveX))
	X_axis = C.create_line(0, MoveY, RozmiarX, MoveY, fill = 'red')
	Y_axis = C.create_line(MoveX, 0, MoveX, RozmiarY, fill = 'red')


Button1 = Tkinter.Button(frame, text ="RIGHT", command = Action1)
Button1.pack(side = LEFT)

Button2 = Tkinter.Button(frame, text ="LEFT", command = Action2)
Button2.pack(side = LEFT)

Button3 = Tkinter.Button(frame, text ="UP", command = Action3)
Button3.pack(side = LEFT)

Button4 = Tkinter.Button(frame, text ="DOWN", command = Action4)
Button4.pack(side = LEFT)

Button5 = Tkinter.Button(frame, text ="IN", command = Action5)
Button5.pack(side = LEFT)

Button6 = Tkinter.Button(frame, text ="OUT", command = Action6)
Button6.pack(side = LEFT)

Button7 = Tkinter.Button(frame, text ="N++", command = Action7)
Button7.pack(side = LEFT)

Button8 = Tkinter.Button(frame, text ="N--", command = Action8)
Button8.pack(side = LEFT)

Button9 = Tkinter.Button(frame, text ="PRINT", command = Action9)
Button9.pack(side = LEFT)

Button10 = Tkinter.Button(frame, text ="INIT", command = Action10)
Button10.pack(side = LEFT)

Label1 = Label(frame, textvariable=var1, relief=RAISED)
var1.set("N = 0")
Label1.pack(side = LEFT)

Label2 = Label(frame, textvariable=var2, relief=RAISED)
var2.set("Zoom = 0")
Label2.pack(side = LEFT)

initialize()

root.mainloop()
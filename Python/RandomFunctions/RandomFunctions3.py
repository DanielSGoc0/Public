from tkinter import *
import tkinter
import copy
import random
import math
import numpy


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
coeffs = []


def new_random_function():
	global N, coeffs

	coeffs = numpy.random.normal(0, 1, N)

def f(x):
	global N, coeffs
	y = 0
	for i in range(N - 1, -1, -1):
		y = y * x
		y = y + coeffs[i]
	return y
	


def evaluate(x):
	global Movex, MoveY

	if x % 2 == 0:
		return x
	else:
		x -= MoveX
		x *= 2**(-Zoom)

		y = f(math.tanh(x)) * 1.0/math.cosh(x)

		return MoveY - y * 2**Zoom



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
	new_random_function()
	paint()
	
def Action9():
	global N, var1
	N = N * 2
	var1.set("N = " + str(N))
	paint()

def Action10():
	global N, var1
	if N > 1:
		N = N//2
		var1.set("N = " + str(N))
	paint()



def paint():
	global MoveX, MoveY, line_new, cord, X_axis, Y_axis, RozmiarX, RozmiarY
	C.delete(ALL)

	for w in range(0, 160):
		cord = [evaluate(i) for i in range(w*12, (w + 1)*12 + 2)]
		line_new = C.create_line(cord)
	#print(MoveY - funkcja(11.5 + 0.3333333 + MoveX))
	X_axis = C.create_line(0, MoveY, RozmiarX, MoveY, fill = 'red')
	Y_axis = C.create_line(MoveX, 0, MoveX, RozmiarY, fill = 'red')

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




#!/usr/bin/python3

from Tkinter import *
import Tkinter
from bigfloat import *

RozmiarX = 1366
RozmiarY = 670

MoveX = RozmiarX/2
MoveY = RozmiarY/2

Zoom = 0
N = 0

root = Tk()
frame = Frame(root)
frame.pack()

var1 = StringVar()
var2 = StringVar()
cord = [0, 0, 0, 0]

seq = []

C = Tkinter.Canvas(root, height=RozmiarY, width=RozmiarX, bg = 'white')
C.pack()

line_old = C.create_line(cord, fill = 'white')
line_new = C.create_line(cord)

X_axis = C.create_line(cord, fill = 'red')
Y_axis = C.create_line(cord, fill = 'red')

def F(x):
	return zeta(2*x)

def intitialize():
	global seq
	with precision(200):
		seq = [BigFloat(F(i)) for i in range(0, 32)]


def Binomial(a, b):
	if b > a:
		return 0
	if b == 0:
		return 1
	if b > a/2:
		return Binomial(a, a - b)
	w = 1
	for i in range(1, b + 1):
		w *= a + 1 - i
		w /= i
	return w

def funkcja(x):
	global Movex, MoveY, N, seq
	if x % 2 == 0:
		return x
	else:
		if len(seq) < 20:
			intitialize()
		x -= MoveX
		x *= 2**(-Zoom)
		x -= 0.666666666666666
		with precision(200):
			# MAIN PART
			# argument to obecne x
			y = BigFloat(0)
			for k in range(0, N + 1):
				z = mul(mul(seq[k + 1], div(BigFloat(N + 1 - x), BigFloat(k + 1 - x))),
					mul(BigFloat(Binomial(N + k + 1, N)), BigFloat(Binomial(N, k))))
				if (N - k) % 2 == 0:
					y = add(y, z)
				else:
					y = sub(y, z)
			for i in range(1, N + 1):
				y = mul(y, div(BigFloat(x - i), BigFloat(x + i)))
			# MAIN PART
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
	global N, var1
	if N < 30:
		N += 1
		var1.set("N = " + str(N))
	paint()

def Action8():
	global N, var1
	if N > 0:
		N -= 1
		var1.set("N = " + str(N))
	paint()

def Action9():
	global N, var1
	if N < 21:
		N += 10
		var1.set("N = " + str(N))
	paint()

def Action10():
	global N, var1
	if N > 0:
		N -= 10
		var1.set("N = " + str(N))
	paint()

def paint():
	global MoveX, MoveY, line_new, cord, X_axis, Y_axis
	C.delete(ALL)
	for w in range(0, 114):
		cord = [funkcja(i) for i in range(w*12, (w + 1)*12 + 2)]
		line_new = C.create_line(cord)
	X_axis = C.create_line(0, MoveY, 1366, MoveY, fill = 'red')
	Y_axis = C.create_line(MoveX, 0, MoveX, 670, fill = 'red')


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

Button9 = Tkinter.Button(frame, text ="N+10", command = Action9)
Button9.pack(side = LEFT)

Button10 = Tkinter.Button(frame, text ="N-10", command = Action10)
Button10.pack(side = LEFT)

Label1 = Label(frame, textvariable=var1, relief=RAISED)
var1.set("N = 0")
Label1.pack(side = LEFT)

Label2 = Label(frame, textvariable=var2, relief=RAISED)
var2.set("Zoom = 0")
Label2.pack(side = LEFT)

root.mainloop()
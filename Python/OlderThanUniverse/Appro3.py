#!/usr/bin/python3

from tkinter import *
import tkinter
from math import *
from scipy.special import *
import random
#from bigfloat import *

euler_gamma = 0.57721566490153286060651209008240243104215933593992
RozmiarX = 1800
RozmiarY = 1000

MoveX = RozmiarX/2
MoveY = RozmiarY/2

Zoom = 0
N = 0
S = 1.0/2
INC = 0

root = Tk()
frame = Frame(root)
frame.pack()

var1 = StringVar()
var2 = StringVar()
var3 = StringVar()
cord = [0, 0, 0, 0]

seq = []

C = tkinter.Canvas(root, height=RozmiarY, width=RozmiarX, bg = 'white')
C.pack()

line_old = C.create_line(cord, fill = 'white')
line_new = C.create_line(cord)

X_axis = C.create_line(cord, fill = 'red')
Y_axis = C.create_line(cord, fill = 'red')

DIFF = False
DOTS = True

def F(x):
	# return random.uniform(0.99, 1.0) / 1.05**x
	# s = 2*x + 2
	# s = 2*x
	s = x + 0.0001
	# return zeta(s)*(s-1)
	return 1/zeta(2*s)
	a = 1
	# return (zeta(s)*a**(2*s) - a**2/(s-1))
	# return a**(2*x)*zeta(2*x) - a/(2*x - 1)
	# return zeta(s)*(1 - 3**(1-s))

def intitialize():
	global seq
	seq = [F(i) for i in range(0, 32)]

def Legendre(x, n):
	R = 0.0
	for i in range(0, n + 1):
		y = (-x) ** i
		y *= Binomial(n + i, i)
		y *= Binomial(n, i)
		R += y
	return R

def Binomial(a, b):
	if b > a:
		return 0
	if b == 0:
		return 1
	if b > a/2:
		return Binomial(a, a - b)
	w = 1
	for i in range(1, b + 1):
		w *= (a + 1 - i)
		w /= i
	return w

def funkcja(x):
	global Movex, MoveY, N, seq, S
	if x % 2 == 0:
		return x
	else:
		#if len(seq) < 20:
		#	intitialize()
		x -= MoveX
		x -= 0.3333333
		x *= 2**(-Zoom)
		# MAIN PART
		# argument to obecne x
		#x = x + 14.134725j/4.0
		y = 0.0
		for k in range(0, N + 1):
			z = seq[k]
			if DIFF:
				z -= F(x)
			z *= (-1) ** k
			z *= Binomial(N + k, N)
			z *= Binomial(N, k)
			z *= x/(k - x)
			y += z
		if DIFF:
			y *= 1
		if not DIFF:
			for i in range(0, N + 1):
				y *= (i - x)/(i + x)
		return MoveY - y * 2**(Zoom + INC)

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
	global Zoom
	Zoom += 1
	var2.set("Zoom = " + str(Zoom))
	paint()

def Action6():
	global Zoom
	Zoom -= 1
	var2.set("Zoom = " + str(Zoom))
	paint()

def Action7():
	global N, var
	if N < 30:
		N += 1
		var1.set("N = " + str(N))
	paint()

def Action8():
	global N, var
	if N > 0:
		N -= 1
		var1.set("N = " + str(N))
	paint()

def Action9():
	global DIFF
	if DIFF:
		DIFF = False
	else:
		DIFF = True
	paint()

def Action10():
	global DOTS
	if DOTS:
		DOTS = False
	else:
		DOTS = True
	paint()

def Action11():
	global INC
	if INC < 100:
		INC += 1
		var3.set("INC = " + str(INC))
	paint()

def Action12():
	global INC
	if INC > -100:
		INC -= 1
		var3.set("INC = " + str(INC))
	paint()

def paint():
	global MoveX, MoveY, line_new, cord, X_axis, Y_axis, RozmiarX, RozmiarY
	C.delete(ALL)
	intitialize()
	for w in range(0, 114):
		cord = [funkcja(i) for i in range(w*12, (w + 1)*12 + 2)]
		line_new = C.create_line(cord)
	#print(MoveY - funkcja(11.5 + 0.3333333 + MoveX))
	X_axis = C.create_line(0, MoveY, RozmiarX, MoveY, fill = 'red')
	Y_axis = C.create_line(MoveX, 0, MoveX, RozmiarY, fill = 'red')
	if DOTS:
		for w in range(0, N + 1):
			PX = w * 2**(Zoom) + MoveX
			PY = MoveY - seq[w] * 2**Zoom
			C.create_oval(PX - 2, PY - 2, PX + 2, PY + 2)


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

Button7 = tkinter.Button(frame, text ="N++", command = Action7)
Button7.pack(side = LEFT)

Button8 = tkinter.Button(frame, text ="N--", command = Action8)
Button8.pack(side = LEFT)

Button9 = tkinter.Button(frame, text ="CHANGE!!!!!!!!", command = Action9)
Button9.pack(side = LEFT)

Button10 = tkinter.Button(frame, text ="Dots", command = Action10)
Button10.pack(side = LEFT)

Button11 = tkinter.Button(frame, text ="INC++", command = Action11)
Button11.pack(side = LEFT)

Button12 = tkinter.Button(frame, text ="INC--", command = Action12)
Button12.pack(side = LEFT)

Label1 = Label(frame, textvariable=var1, relief=RAISED)
var1.set("N = 0")
Label1.pack(side = LEFT)

Label2 = Label(frame, textvariable=var2, relief=RAISED)
var2.set("Zoom = 0")
Label2.pack(side = LEFT)

Label3 = Label(frame, textvariable=var3, relief=RAISED)
var3.set("INC = 0")
Label3.pack(side = LEFT)

root.mainloop()

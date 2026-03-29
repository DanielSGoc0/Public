# This program is meant for checking the hypothesis that e^(-ax^2) is optimal

from tkinter import *
import tkinter
import copy
import random
import math
import numpy
import scipy
# import mpmath


RozmiarX = 1800
RozmiarY = 1000

MoveX = RozmiarX/2
MoveY = RozmiarY/2
Zoom = 8

root = Tk()
frame = Frame(root)
frame.pack()

var1 = StringVar()
var2 = StringVar()
var3 = StringVar()
var4 = StringVar()

C = tkinter.Canvas(root, height=RozmiarY, width=RozmiarX, bg = 'white')
C.pack()
NODES = []
L = 32
N = 16
RES = 5

def sign(x):
	if x > 0:
		return 1
	else:
		return 0
	
def sech(x):
	if abs(x) > 100:
		return 0
	else:
		return 1.0/numpy.cosh(x)

def f(x, a):
	global L, RES
	# phi1 = lambda t: 1
	# phi1 = lambda t: 1.0/(1 + t*t)
	# phi1 = lambda t: math.exp(-abs(t))
	# phi1 = lambda t: t**2 * math.exp(-t*t/2)
	# phi1 = lambda t: sech(t) * math.exp(-t*t)
	# phi1 = lambda t: max(0, 1 - abs(2*t)) 
	# phi1 = lambda t: max(0, 1 - abs(t))
	# phi1 = lambda t: math.exp(-t**4)
	phi1 = lambda t: t**2 * math.exp(-t*t/2)
	# phi1 = lambda t: math.exp(-t*t/2)
	# phi2 = lambda t: sign(1 - abs(t))
	# phi2 = lambda t: 1.0/(1 + t*t)**20
	# phi2 = lambda t: sign(1 - abs(t - 2)) + sign(1 - abs(t + 2))
	# phi2 = lambda t: 1.0/2*(math.exp(-t*t/10) + math.exp(-t*t*10))
	# phi2 = lambda t: numpy.sinc(t)**2
	phi2 = lambda t: math.exp(-t*t/2)
	# phi2 = lambda t: 1.0/(1 + t*t)
	# phi2 = lambda t: math.exp(-abs(t))
	# phi2 = lambda t: 1

	# a = L / 2**RES
	# print(a)
	res = scipy.integrate.quad(lambda t: phi1(t)**a * phi2(t)**(1 - a) * math.cos(4*x*t), -numpy.inf, numpy.inf)
	return res[0]
	# if abs(res[0]) < 0.00000000001:
	# 	return -10
	# elif res[0] < 0:
	# 	return 10
	# else:
	# 	return math.log(res[0])


def new_plot():
	global N, NODES
	NODES.clear()

	for k in range(-3*N, 3*N + 1):
		NODES.append((k/N, f(k/N, L / 2**RES)))

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
	global L, var3
	if L < 2**RES:
		L += 1
	var3.set("L = " + str(L))
	new_plot()
	paint()

def Action8():
	global L, var3
	if L > 0:
		L -= 1
	var3.set("L = " + str(L))
	new_plot()
	paint()
	
def Action9():
	global N, var1
	N = N * 2
	var1.set("N = " + str(N))
	new_plot()
	paint()

def Action10():
	global N, var1
	if N > 1:
		N = N // 2
		var1.set("N = " + str(N))
	new_plot()
	paint()

def Action11():
	global RES, L, var4
	RES += 1
	L *= 2
	var3.set("L = " + str(L))
	var4.set("RES = " + str(RES))
	new_plot()
	paint()

def Action12():
	global RES, L, var4
	if RES > 0:
		RES -= 1
		L = L//2
		var3.set("L = " + str(L))
		var4.set("RES = " + str(RES))
	new_plot()
	paint()



def paint():
	global MoveX, MoveY, X_axis, Y_axis, RozmiarX, RozmiarY, NODES
	C.delete(ALL)

	X_axis = C.create_line(0, MoveY, RozmiarX, MoveY, fill = 'red')
	Y_axis = C.create_line(MoveX, 0, MoveX, RozmiarY, fill = 'red')

	for i in range(len(NODES) - 1):
		PX0 = NODES[i][0] * 2**(Zoom) + MoveX
		PY0 = MoveY - NODES[i][1] * 2**Zoom
		PX1 = NODES[i + 1][0] * 2**(Zoom) + MoveX
		PY1 = MoveY - NODES[i + 1][1] * 2**Zoom

		C.create_line(PX0, PY0, PX1, PY1, fill = 'black')

	# global m, M, I
	# C.create_oval(m * 2**(Zoom) + MoveX - 2, MoveY - 2, m * 2**(Zoom) + MoveX + 2, MoveY + 2)
	# C.create_oval(I * 2**(Zoom) + MoveX - 2, MoveY - 2, I * 2**(Zoom) + MoveX + 2, MoveY + 2)
	# C.create_oval(M * 2**(Zoom) + MoveX - 2, MoveY - 2, M * 2**(Zoom) + MoveX + 2, MoveY + 2)



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

Button7 = tkinter.Button(frame, text ="L++", command = Action7)
Button7.pack(side = LEFT)

Button8 = tkinter.Button(frame, text ="L--", command = Action8)
Button8.pack(side = LEFT)

Button9 = tkinter.Button(frame, text ="N++", command = Action9)
Button9.pack(side = LEFT)

Button10 = tkinter.Button(frame, text ="N--", command = Action10)
Button10.pack(side = LEFT)

Button11 = tkinter.Button(frame, text ="RES++", command = Action11)
Button11.pack(side = LEFT)

Button12 = tkinter.Button(frame, text ="RES--", command = Action12)
Button12.pack(side = LEFT)



Label1 = Label(frame, textvariable=var1, relief=RAISED)
var1.set("N = 16")
Label1.pack(side = LEFT)

Label2 = Label(frame, textvariable=var2, relief=RAISED)
var2.set("Zoom = 0")
Label2.pack(side = LEFT)

Label3 = Label(frame, textvariable=var3, relief=RAISED)
var3.set("L = 16")
Label3.pack(side = LEFT)

Label4 = Label(frame, textvariable=var4, relief=RAISED)
var4.set("RES = 5")
Label4.pack(side = LEFT)

new_plot()
paint()
root.mainloop()


# for a in range(2, 3):
# 	print("a = ", a)
# 	for i in range(1, 21):
# 		print(scipy.integrate.quad(lambda t: f(t, 1) * math.exp(a*t/3.0) * math.sin(t/3.0), -i, i)[0])
		# print(scipy.integrate.quad(lambda t: math.exp(a*t/3.0 - t*t) * math.cos(a*t/3.0), -i, i)[0])


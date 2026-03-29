# Original plotter. Don't use!

from tkinter import *
import tkinter
import copy
import random
import math
import numpy as np
from numpy import *
import scipy
import scipy.integrate
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

C = tkinter.Canvas(root, height=RozmiarY, width=RozmiarX, bg = 'white')
C.pack()
N = 1024
NODES = []

res = []
a = 0.3
c = 0.7
I1 = 0.3
I2 = 0.7
gamma = 1
A = 0.5

def magic():
	global a, c, I1, I2, gamma

	tab = [random.random() for _ in range(4)]

	tab = sorted(tab)

	a = tab[0]
	c = tab[3]
	I1 = tab[1]
	I2 = tab[2]
	gamma = random.random()
	# M = m + 0.000001



def F(x):
	global A
	# return math.log(x)
	# return -x*math.log(x)
	# return -x*math.log(x) - (1 - x)*math.log(1 - x)
	return x**A + (1 - x)**A
	# return x**a
	# return -x**a
	# return - x**a - (1 - x)**a
	# if x < 1.0/3:
	# 	return (x - 1.0/3) * 10
	# elif x < 2.0/3:
	# 	return 0
	# else:
	# 	return -(x - 2.0/3) * 1

def dF(x):
	global A
	# return 1/x
	# return -math.log(x) - 1
	# return math.log((1 - x)/x)
	return A*x**(A - 1) - A*(1-x)**(A - 1)
	# return a*x**(a - 1)
	# return -a*x**(a - 1)
	# return -a * x**(a - 1) + a *(1-x)**(a - 1)
	# if x < 1.0/3:
	# 	return 10
	# elif x < 2.0/3:
	# 	return 0
	# else:
	# 	return -1

def ddF(x):
	global A
	# return -1/(x*x)
	# return -1/x
	# return -1/(x*(1-x))
	return A*(A - 1)*x**(A - 2) + A*(A - 1)*(1-x)**(A - 2)


def optimal_m(v):
	a = 0.0
	b = 1.0
	for _ in range(70):
		h = (a + b)/2
		if dF(h) > v:
			a = h
		else:
			b = h
	return (a + b)/2


def optimum(a, b):
	S = (F(b) - F(a))/(b - a)
	m = optimal_m(S)
	# return m
	return F(m) - ((b - m) * F(a) + (m - a) * F(b))/(b - a)

def difference(a, b, c):
	return optimum(a, c) - (optimum(a, b) + optimum(b, c))


def H(x):
	return -x*math.log(x) - (1 - x)*math.log(1-x)
	# return -x*math.log(x)

def I(a, b):
	return 1/(1 + math.exp((H(b) - H(a))/(b - a)))

def w(x):
	return math.log(1 + x * math.exp(math.log(1 - x) * (1 - x) / x))

def mean(x, y, z):
	return H(y) - ((y - x) * H(z) + (z - y) * H(x))/(z - x)

def second_derivative(a, b):
	I = optimal_m((F(b) - F(a))/(b - a))
	R = (dF(b) - dF(I))/(b - a)
	T = (I - a)/(b - a)

	return -(R*R/ddF(I) - 2*R*T + T*ddF(b))
	# return R*R/ddF(I) + T*ddF(b)

def error(X, S, STEPS):
	X = np.array(X)
	N = X.shape[0]
	h = 1.0 * S/STEPS

	M = np.zeros((N, N), dtype='double')
	for i in range(N):
		for j in range(N):
			M[i][j] = math.exp(-np.dot(X[i] - X[j], X[i] - X[j]) / 2.0)
	V = np.linalg.cholesky(M)

	# just run a forward Euler
	g = np.zeros((N), dtype='double')
	w = np.zeros((N), dtype='double')

	for _ in range(STEPS):
		for i in range(N):
			g[i] += h * math.exp(np.dot(X[i], w - X[i]/2.0))
		w = scipy.linalg.solve_triangular(V, g, lower=True)
	
	return -2*math.log(1 - S*S/2) - np.dot(w, w)

def f(x):
	global a, c, I1, I2, gamma
	# d = c - a

	# if x < 0.01 or x > 0.99:
	# 	return 0
	
	# a = 0.0001
	# c = 0.9999

	# b = x
	# h = 0.0001
	# I3 = gamma*I2 + (1-gamma)*I1

	# a = 0.00001
	# if x < 0.0001 or b + x > 0.9999:
	# 	return 0

	# return optimum(a, b)
	# return (optimum(a, b - h) + optimum(a, b + h) - 2*optimum(a, b))/(h*h) / 10
	# return second_derivative(a, b)/10
	# I = optimal_m((F(b) - F(a))/(b - a))
	# return (optimum(x, x - h/2) + optimum(x, x + 3*h/2)-2*optimum(x, x + h/2))/(h*h)
	# return (optimum(a, b + h) - optimum(a, b))/h
	# return (I - a)/(b - a)*(dF(I) - dF(b))
	# return (optimum(a, b + h) - optimum(a, b))/h - (I - a)/(b - a)*(dF(I) - dF(b))
	# return (optimum(a, b - h) + optimum(a, b + h)-2*optimum(a, b))/(h*h) - second_derivative(a, b)
	# return (mean(a - a + x, I3 - a + x, c - a + x) - (mean(a - a + x, I1 - a + x, b - a + x) + mean(b - a + x, I2 - a + x, c - a + x)))
	# I2 = 1.0/(1 + math.exp((H(a + x) - H(a))/x))
	# return (I2 - a)/x
	# Y = 1/(b * (1 + math.exp(H(b)/b)))
	# return x * math.log(1/(b + x*Y) - 1)
	# return x * math.log(1/(b + x*Y) - 1) - H(b + x) + H(b)
	# A = Y
	# B = b
	# return (A*x)/((A*x + B - 1)*(A*x + B)) - math.log((A*x + B)/(x + B)) - ((1 - B - x)/(1 - A*x - B) - 1)
	# return (A*x)/((A*x + B - 1)*(A*x + B)) - math.log((A*x + B)/(x + B)) - ((1 - B - x)/(1 - A*x - B) - 1)
	# return -A/((1 - A*x - B)*(A*x + B)) + (1 - A)/(x + B) * (A*x + x + 2*B)/(A*x + B) * 1.0/2 + (1 - A)/(1 - A*x - B)
	# return -A + (1 - A)/(x + B) * (A*x + x + 2*B)*(1 - A*x - B) * 1.0/2 + (1 - A) * (A*x + B)
	# return -2*A*(x + B) + (1 - A)*(x + B)*(1 + A*x + B) + (1 - A)*(A*x + B)*(1 - A*x - B)
	# return (1 + B - math.sqrt(2)* (1 + B))/(B - 1) - A
	# return (I(b, b + x) - b)/x - 1/2.0
	# return 1 + B*(1 - A)**2 - 2*A - A**2

	# func = lambda t: abs(np.sinc(t))**3
	# print(x)
	# val = scipy.integrate.quad(lambda t: t**x * math.exp(func(t)), 0, np.inf)
	# val = scipy.integrate.quad(lambda t: t**x * math.exp(func(t)), 0, 1)

	# phi = scipy.integrate.quad(lambda t: func(t) * math.cos(4*x*t), -np.inf, np.inf)
	# M = scipy.integrate.quad(lambda t: math.exp(func(t) + 4*x*t) + math.exp(func(t) - 4*x*t), 0, np.inf)
	# print(res)
	# return res
	# return phi[0] * M[0]
	# return math.log(val[0] / scipy.special.gamma(x + 1))
	# return math.log(val[0])
	# return phi[0]

	# return 2 * scipy.special.Gamma(7.0/6) HypergeometricPFQ[{}, {1/3, 1/2, 2/3, 5/6}, -t^6/46656] - (t^2 (72 Sqrt[Pi] HypergeometricPFQ[{}, {2/3, 5/6, 7/6, 4/3}, -t^6/46656] + t^2 Gamma[-1/6] HypergeometricPFQ[{}, {7/6, 4/3, 3/2, 5/3}, -t^6/46656]))/432

	# q = 5.0
	# VAL1 = -math.log(1 + q**x) + 1.0/(1 + q**x) * (q**x * math.log(q) * (1-q)/(1 + q) + (q + q**x)/(1 + q) * math.log(1 + q**(x-1)) + (1 + q**(x + 1))/(1 + q) * math.log(1 + q**(x +1)))
	# return 10*VAL
	# VAL2 = 1/(1 + q**x)*(-(1 + q**x) * H(1/(1 + q**x)) + (q + q**x)/(1 + q) * H(1/(1 + q**(x-1))) + (1 + q**(x+1))/(1 + q) * H(1/(1 + q**(x + 1)))) 
	# return 10*VAL2 
	# VAL3 = -math.log(1 + x) + 1.0/(1 + x) * (x * math.log(q) * (1-q)/(1 + q) + (q + x)/(1 + q) * math.log(1 + x/q) + (1 + x*q)/(1 + q) * math.log(1 + x*q))
	# return VAL3

	# return error([[0.0, 0.0, 0.0], [x + 0.001, 0.0, 0.0], [0.6563 + x, 1.0 - 3*x, 0.0]], 1, 1000)

	# v = scipy.integrate.quad(lambda t: math.log(1 + math.exp(-t*t)) * math.cos(5*x*t), 0, np.inf)[0]
	# print(5*x, v)


	# h = np.exp(-3*x)
	# z0 = 1.072757746848
	# z1 = 0.479279265336
	# z2 = 2.048374923749
	# x0 = 0.0
	# x1 = z0 * h - z0*z0*z0*h*h*h
	# x2 = z0 * h + z0 * z1 * h*h + z0*(z1*z1 - z0*z0/2) * h*h*h + z0*z1*(z1**2 - 7.0/4 * z0**2) * h**4
	# # x2 = z0 * h + z0 * z1 * h*h + z0*(z1*z1 - z0*z0/2) * h*h*h
	# X = np.array([[x0], [x1], [x2]])
	# S = np.zeros((4, 4))
	# for i in range(3):
	# 	for j in range(3):
	# 		S[i][j] = np.exp(-np.dot(X[i] - X[j], X[i] - X[j]) / 2.0)
	# C = np.zeros((4, 4))
	# C[0][0] = 1.0
	# C[1][0] = S[1][0]
	# C[1][1] = np.sqrt(1.0 - S[1][0]**2) * np.sign(X[1][0])
	# C[2][0] = S[2][0]
	# C[2][1] = (S[2][1] - C[2][0] * C[1][0])/C[1][1]
	# C[2][2] = np.sqrt(1 - C[2][1]**2 - C[2][0]**2)
	# # C[3][0] = S[3][0]
	# # C[3][1] = (S[3][1] - C[3][0] * C[1][0])/C[1][1]
	# # C[3][2] = (S[3][2] - C[3][0] * C[2][0] - C[3][1] * C[2][1]) / C[2][2]

	# print((x2 - h*(C[2][0] * z0 + C[2][1] * z1)) / h**5)
	# return (x2 - h*(C[2][0] * z0 + C[2][1] * z1)) / h**5
	# # print((x2 - h*(C[2][0] * z0 + C[2][1] * z1 + C[2][2] * z2))**2 / h**8)
	# # return (x2 - h*(C[2][0] * z0 + C[2][1] * z1 + C[2][2] * z2))**2 / h**8

	return 2*floor(3*x) - 3*floor(2*x)


def generate_list():
	global N, res
	res = []
	for k in range(3*N + 1):
		func = lambda t: 100*t
		# val = scipy.integrate.quad(lambda t: t**(k/N) * math.exp(func(t)), 0, np.inf)
		val = scipy.integrate.quad(lambda t: t**(k/N) * math.exp(func(t)), 0, 0.2)
		# res.append(val[0] / scipy.special.gamma((k/N + 1)/2))
		res.append(val[0])


def new_plot():
	global N, NODES, res
	NODES.clear()

	# generate_list()
	for k in range(-2*N, 2*N + 1):
		NODES.append((k/N, f(k/N)))
	# for k in range(2, 3*N):
	# 	NODES.append((k/N, N*N*math.log(res[k - 1] * res[k + 1] / res[k]**2)))
	# print(res[0] * res[2*N] / res[N]**2)
	# print(res[N] * res[3*N] / res[2*N]**2)

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
	new_plot()
	paint()

def Action8():
	magic()
	paint()
	
def Action9():
	global N, var1
	N = N * 2
	var1.set("N = " + str(N))
	paint()

def Action10():
	global N, var1
	if N > 1:
		N = N // 2
		var1.set("N = " + str(N))
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

Button7 = tkinter.Button(frame, text ="NEW PLOT", command = Action7)
Button7.pack(side = LEFT)

Button8 = tkinter.Button(frame, text ="???", command = Action8)
Button8.pack(side = LEFT)

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




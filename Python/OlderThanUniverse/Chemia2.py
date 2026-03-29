
from Tkinter import *
import Tkinter
from math import *
import random
#from bigfloat import *

RozmiarX = 1366
RozmiarY = 670

MoveX = RozmiarX/2
MoveY = RozmiarY/2

TRIES = 200
AMT = 10
N = 10
DATA_SIZE = 13

STALA = 0.01

root = Tk()
frame = Frame(root)
frame.pack()

var1 = StringVar()
var2 = StringVar()
var3 = StringVar()
var4 = StringVar()
var5 = StringVar()
cord = [0, 0, 0, 0]

shown = True

seq = []

px = []
py = []

generator = []
optimum = []
present = []

C = Tkinter.Canvas(root, height=RozmiarY, width=RozmiarX, bg = 'white')
C.pack()

line_old = C.create_line(cord, fill = 'white')
line_new = C.create_line(cord)

X_axis = C.create_line(cord, fill = 'red')
Y_axis = C.create_line(cord, fill = 'red')

#================================================   KONWERSJA  ==========================================================

def convert():
	tab1 = [0 for i in range(0, DATA_SIZE)]
	tab2 = [0 for i in range(0, DATA_SIZE)]
	#tab3 = [0 for i in range(0, DATA_SIZE)]
	#tab4 = [0 for i in range(0, DATA_SIZE)]
	plik = open('Chemia/text2.txt')
	s = plik.read()
	#print s
	global px, py
	number = 1
	ind = 0
	dot = -1
	x = 0
	START = False
	minus = True
	for i in range(0, len(s)):
		if 0 <= ord(s[i]) - 48 and ord(s[i]) - 48 <= 9:
			if dot >= 0:
				dot = dot + 1
			if START:
				x *= 10
				x += ord(s[i]) - 48
			else:
				START = True
				x = ord(s[i]) - 48
		elif s[i] == '-':
			minus = True
		elif s[i] == '.':
			dot = 0
		elif START:
			minus = False
			START = False
			x *= 1.0000000000001
			if dot >= 0:	
				x /= 10**dot
			if minus:
				x = -x
			#print(x)
			if number == 1:
				tab1[ind] = x
			else:
				tab2[ind] = x
				number = 0
				ind = ind + 1
			#elif number == 2:
			#	tab2[ind] = x
			#elif number == 3:
			#	tab3[ind] = x
			#else:
			#	tab4[ind] = x
			#	number = 0
			#	ind = ind + 1
			number = number + 1
			x = 0
			dot = -1
	px = [0 for i in range(0, DATA_SIZE)]
	py = [0 for i in range(0, DATA_SIZE)]
	for i in range(0, DATA_SIZE):
		px[i] = tab1[i]
		py[i] = tab2[i] * 3.0
		#print(str(tab2[i]) + "  " + str(tab3[i]))
	N = DATA_SIZE

#=============================================================================================================================


def value_generator(x):
	global generator
	return generator[0] * (generator[1] / (x / generator[2] + 1) + generator[3])

def distance_generator():
	global px, py
	D = 0.0
	for i in range(0, N):
		dy = value_generator(px[i]) - py[i]
		D += dy*dy
	return D

def gen():
	global N, generator, px, py
	generator = [random.random() for i in range(0, 4)]
	px = [0 for i in range(0, N)]
	py = [0 for i in range(0, N)]
	for i in range(0, N):
		px[i] = random.random()
		py[i] = value_generator(px[i]) + 0.2 * (0.5 - random.random())
	
	print("dist_gen:" + str(distance_generator()))

def value_optimum(x):
	global optimum
	return optimum[0] * (optimum[1] / (x / optimum[2] + 1) + optimum[3])

def distance_optimum():
	global px, py
	D = 0.0
	for i in range(0, N):
		dy = value_optimum(px[i]) - py[i]
		D += dy*dy
	return D

def value_present(x):
	global present
	return present[0] * (present[1] / (x / present[2] + 1) + present[3])

def distance_present():
	global px, py
	D = 0.0
	for i in range(0, N):
		dy = value_present(px[i]) - py[i]
		D += dy*dy
	return D

def da():
	global present, px, py, N
	DA = 0.0
	for i in range(0, N):
		DA += (value_present(px[i]) - py[i]) * (present[1] / (px[i] / present[2] + 1) + present[3])
	return 2*DA

def db():
	global present, px, py, N
	DB = 0.0
	for i in range(0, N):
		DB += (value_present(px[i]) - py[i]) * present[0] / (px[i] / present[2] + 1)
	return 2*DB

def dc():
	global present, px, py, N
	DC = 0.0
	for i in range(0, N):
		DC += (value_present(px[i]) - py[i]) * present[0] * present[1] * px[i] / ((px[i] + present[2]) * (px[i] + present[2]))
	return 2*DC

def dd():
	global present, px, py, N
	DD = 0.0
	for i in range(0, N):
		DD += (value_present(px[i]) - py[i]) * present[0]
	return 2*DD

def find_optimum():
	global optimum, present, TRIES, px, py
	optimum = [random.random() for i in range(0, 4)]
	for i in range(0, TRIES):
		present = [random.random() for j in range(0, 4)]
		#print(present[0])
		present[0] = 0.000663231095446 * (0.5 + present[0])
		present[1] = -1764.12615556 * (0.5 + present[1])
		present[2] = 0.00152432913726 * (0.5 + present[2])
		present[3] = 6395.97734583 * (0.5 + present[3])
		for j in range(0, 2**AMT):
			Da = da()
			Db = db()
			Dc = dc()
			Dd = dd()
			#if j % 10000 == 0:
			#	print("=============================")
			#	print(Da)
			#	print(Db)
			#	print(Dc)
			#	print(Dd)
			#	print("------       " + str(j))
			#	print(present[0])
			#	print(present[1])
			#	print(present[2])
			#	print(present[3])
			#print(str(Da) + "  " + str(Db) + "  " + str(Dc) + "  " + str(Dd))
			#print(str(generator[0]) + "\t" + str(generator[1]) + "\t" + str(generator[2]) + "\t" + str(generator[3]))
			#print(str(present[0]) + "\t" + str(present[1]) + "\t" + str(present[2]) + "\t" + str(present[3]))
			# if fabs(Da) >= fabs(Db) and fabs(Da) >= fabs(Dc) and fabs(Da) >= fabs(Dd):
				#print("Da")
			#	present[0] -= Da / 10.0
			#if fabs(Db) >= fabs(Da) and fabs(Db) >= fabs(Dc) and fabs(Db) >= fabs(Dd):
				#print("Db")
			#	present[1] -= Db / 10.0
			#if fabs(Dc) >= fabs(Da) and fabs(Dc) >= fabs(Db) and fabs(Dc) >= fabs(Dd):
				#print("Dc")
			#	present[2] -= Dc / 10.0
			#if fabs(Dd) >= fabs(Da) and fabs(Dd) >= fabs(Db) and fabs(Dd) >= fabs(Dc):
				#print("Dd")
			#	present[3] -= Dd / 10.0
			if fabs(present[0]) * STALA < fabs(Da):
				if Da > 0:
					present[0] = present[0] * (1 - STALA)
				else:
					present[0] = present[0] * (1 + STALA)
			else:
				present[0] -= Da * STALA
			if fabs(present[1]) * STALA < fabs(Db):
				if Db > 0:
					present[1] = present[1] * (1 - STALA)
				else:
					present[1] = present[1] * (1 + STALA)
			else:
				present[1] -= Db * STALA
			if fabs(present[2]) * STALA < fabs(Dc):
				if Dc > 0:
					present[2] = present[2] * (1 - STALA)
				else:
					present[2] = present[2] * (1 + STALA)
			else:
				present[2] -= Dc * STALA
			if fabs(present[3]) * STALA < fabs(Dd):
				if Dd > 0:
					present[3] = present[3] * (1 - STALA)
				else:
					present[3] = present[3] * (1 + STALA)
			else:
				present[3] -= Dd * STALA
		if distance_present() < distance_optimum():
			optimum = [present[i] for i in range(0, 4)]
	#var4.set("ER_GEN = " + str(distance_generator()))
	var5.set("ER_OPT = " + str(distance_optimum()))


def Action1():
	global N
	if N > 5:
		N -= 1
		var1.set("N = " + str(N))
	paint()

def Action2():
	global N
	N += 1
	var1.set("N = " + str(N))
	paint()

def Action3():
	gen()
	find_optimum()
	paint()

def Action4():
	global shown
	if shown:
		var2.set("optimum")
		shown = False
	else:
		var2.set("generator")
		shown = True
	paint()

def Action5():
	global AMT
	AMT -= 1
	var3.set("AMT = " + str(AMT))
	find_optimum()
	paint()

def Action6():
	global AMT
	AMT += 1
	var3.set("AMT = " + str(AMT))
	find_optimum()
	paint()

def Action7():
	convert()
	find_optimum()
	paint()

def funkcja_gen(x):
	if x % 2 == 0:
		return x
	else:
		x /= 1366.0
		# MAIN PART
		# argument to obecne x
		y = 100*value_generator(x)
		# MAIN PART
		return 670 - y

def funkcja_opt(x):
	if x % 2 == 0:
		return x
	else:
		x /= 1366.0
		# MAIN PART
		# argument to obecne x
		y = 100*value_optimum(x)
		# MAIN PART
		return 670 - y

def paint():
	global MoveX, MoveY, line_new, cord, X_axis, Y_axis, optimum, generator
	C.delete(ALL)
	if shown:
		if(len(generator) < 2):
			generator = [random.random() for i in range(0, 4)]
		for w in range(0, 114):
			cord = [funkcja_gen(i) for i in range(w*12, (w + 1)*12 + 2)]
			line_new = C.create_line(cord)
	else:
		if(len(optimum) < 2):
			optimum = [random.random() for i in range(0, 4)]
		#print(optimum[0])
		for w in range(0, 114):
			cord = [funkcja_opt(i) for i in range(w*12, (w + 1)*12 + 2)]
			line_new = C.create_line(cord)
	for i in range(0, len(px)):
		#print(str(px[i]) + "  " + str(py[i]))
		C.create_oval(px[i]*1366, 670 - py[i]*100, px[i]*1366, 670 - py[i]*100, width = 3, fill = 'red')
	X_axis = C.create_line(0, 670, 1366, 670, fill = 'red')
	Y_axis = C.create_line(2, 0, 2, 670, fill = 'red')
	print("=============  WYNIK  ==================")
	print(optimum[0])
	print(optimum[1])
	print(optimum[2])
	print(optimum[3])


Button1 = Tkinter.Button(frame, text ="N--", command = Action1)
Button1.pack(side = LEFT)

Button2 = Tkinter.Button(frame, text ="N++", command = Action2)
Button2.pack(side = LEFT)

Button5 = Tkinter.Button(frame, text ="AMT--", command = Action5)
Button5.pack(side = LEFT)

Button6 = Tkinter.Button(frame, text ="AMT++", command = Action6)
Button6.pack(side = LEFT)

Button3 = Tkinter.Button(frame, text ="Initialize", command = Action3)
Button3.pack(side = LEFT)

Button4 = Tkinter.Button(frame, text ="swap", command = Action4)
Button4.pack(side = LEFT)

Button7 = Tkinter.Button(frame, text ="data", command = Action7)
Button7.pack(side = LEFT)

Label1 = Label(frame, textvariable=var1, relief=RAISED)
var1.set("N = 10")
Label1.pack(side = LEFT)

Label3 = Label(frame, textvariable=var3, relief=RAISED)
var3.set("AMT = 10")
Label3.pack(side = LEFT)

Label2 = Label(frame, textvariable=var2, relief=RAISED)
var2.set("generator")
Label2.pack(side = LEFT)

Label4 = Label(frame, textvariable=var4, relief=RAISED)
var4.set("ER_GEN = ?")
Label4.pack(side = LEFT)

Label5 = Label(frame, textvariable=var5, relief=RAISED)
var5.set("ER_OPT = ?")
Label5.pack(side = LEFT)

root.mainloop()
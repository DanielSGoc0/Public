
from Tkinter import *
import Tkinter
from math import *
import random

RozmiarX = 1366
RozmiarY = 670

MoveX = RozmiarX/2
MoveY = RozmiarY/2

TRIES = 100
AMT = 10
N = 13
DATA_SIZE = 13
SSS = 0.05
TTT = 0.5

LC0 = False

root = Tk()
frame = Frame(root)
frame.pack()

var1 = StringVar()
var2 = StringVar()
var3 = StringVar()
var4 = StringVar()
var5 = StringVar()
cord = [0, 0, 0, 0]

TRYB = False

seq = []

px = []
py = []

generator = []
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
		py[i] = tab2[i]
		#print(str(tab2[i]) + "  " + str(tab3[i]))
	N = DATA_SIZE

#=============================================================================================================================

def value_present1(k, x):
	global present
	return (present[k][0] * x * x + present[k][1])/(x*x + present[k][2] * x + present[k][3])

def distance_present1(k):
	global px, py
	D = 0.0
	for i in range(0, N):
		dy = value_present1(k, px[i]) - py[i]
		D += dy*dy
	return D

def da1(k):
	global present, px, py, N
	DA = 0.0
	for i in range(0, N):
		DA += px[i]*px[i]/(px[i]*px[i] + present[k][2] * px[i] + present[k][3]) * (value_present1(k, px[i]) - py[i])
	return 2*DA

def db1(k):
	global present, px, py, N
	DB = 0.0
	for i in range(0, N):
		DB += 1.0/(px[i]*px[i] + present[k][2] * px[i] + present[k][3]) * (value_present1(k, px[i]) - py[i])
	return 2*DB

def dc1(k):
	global present, px, py, N
	DC = 0.0
	for i in range(0, N):
		DC -= px[i] * (present[k][0] * px[i]*px[i] + present[k][1])/((px[i]*px[i] + present[k][2] * px[i] + present[k][3])**2) * (value_present1(k, px[i]) - py[i])
	return 2*DC

def dd1(k):
	global present, px, py, N
	DD = 0.0
	for i in range(0, N):
		DD -= (present[k][0] * px[i]*px[i] + present[k][1])/((px[i]*px[i] + present[k][2] * px[i] + present[k][3])**2) * (value_present1(k, px[i]) - py[i])
	return 2*DD

def value_present2(k, x):
	global present
	return present[k][0] * (present[k][1] / (x / present[k][2] + 1) + present[k][3])

def distance_present2(k):
	global px, py
	D = 0.0
	for i in range(0, N):
		dy = value_present2(k, px[i]) - py[i]
		D += dy*dy
	return D

def da2(k):
	global present, px, py, N
	DA = 0.0
	for i in range(0, N):
		DA += (value_present2(k, px[i]) - py[i]) * (present[k][1] / (px[i] / present[k][2] + 1) + present[k][3])
	return 2*DA

def db2(k):
	global present, px, py, N
	DB = 0.0
	for i in range(0, N):
		DB += (value_present2(k, px[i]) - py[i]) * present[k][0] / (px[i] / present[k][2] + 1)
	return 2*DB

def dc2(k):
	global present, px, py, N
	DC = 0.0
	for i in range(0, N):
		DC += (value_present2(k, px[i]) - py[i]) * present[k][0] * present[k][1] * px[i] / ((px[i] + present[k][2]) * (px[i] + present[k][2]))
	return 2*DC

def dd2(k):
	global present, px, py, N
	DD = 0.0
	for i in range(0, N):
		DD += (value_present2(k, px[i]) - py[i]) * present[k][0]
	return 2*DD

def sortLast(val): 
    return val[4]

# ==============================================================  OPTIMUM  ===============================================================================

def find_optimum():
	global present, TRIES, px, py
	present = [[random.random() for j in range(0, 5)] for i in range(0, TRIES)]
	for j in range(0, TRIES):
		if LC0:
			present[j][0] = 1.000000000001 * random.random() / 10**(3 + 2*random.random())
		else:
			present[j][0] = 0.0000498
		present[j][2] = 1.000000000001 * random.random() / 10**(2 + 2*random.random())
		present[j][3] = 1.000000000001 * 10**(5.5 - 2*random.random())
		present[j][1] = present[j][3] * (1 - 2*random.random())
		#if LC0:
		#	present[j][0] = 10**(-5) * 10**(2*random.random() - 1)
		#else:
		#	present[j][0] = 0.00004
		#present[j][1] = 1.00000000000001 * 10**3 * (5 - 10*random.random())
		#present[j][2] = 1.00000000000001 * 10**(-2) * 10**(2*random.random() - 1)
		#present[j][3] = 1.00000000000001 * 10**3 * 5 * random.random()	
		#print("======================================")
		#print(present[j][0])
		#print(present[j][1])
		#print(present[j][2])
		#print(present[j][3])
	for j in range(0, 2**AMT):
		STALA = SSS/sqrt(j + 10)
		for i in range(0, TRIES):
			# ulepszanie
			Da = 0.001
			Db = 0.001
			Dc = 0.001
			Dd = 0.001
			if TRYB:
				Da = da1(i)
				Db = db1(i)
				Dc = dc1(i)
				Dd = dd1(i)
			else:
				Da = da2(i)
				Db = db2(i)
				Dc = dc2(i)
				Dd = dd2(i)
			if LC0:
				if fabs(present[i][0]) * STALA < fabs(Da):
					if Da > 0:
						present[i][0] = present[i][0] * (1 - 2*random.random() * STALA)
					else:
						present[i][0] = present[i][0] * (1 + 2*random.random() * STALA)
				else:
					present[i][0] -= Da * STALA
			if fabs(present[i][1]) * STALA < fabs(Db):
				if Db > 0:
					present[i][1] = present[i][1] * (1 - 2*random.random() * STALA)
				else:
					present[i][1] = present[i][1] * (1 + 2*random.random() * STALA)
			else:
				present[i][1] -= Db * STALA
			if fabs(present[i][2]) * STALA < fabs(Dc):
				if Dc > 0:
					present[i][2] = present[i][2] * (1 - 2*random.random() * STALA)
				else:
					present[i][2] = present[i][2] * (1 + 2*random.random() * STALA)
			else:
				present[i][2] -= Dc * STALA
			if fabs(present[i][3]) * STALA < fabs(Dd):
				if Dd > 0:
					present[i][3] = present[i][3] * (1 - 2*random.random() * STALA)
				else:
					present[i][3] = present[i][3] * (1 + 2*random.random() * STALA)
			else:
				present[i][3] -= Dd * STALA
			if TRYB:
				present[i][4] = distance_present1(i)
			else:
				present[i][4] = distance_present2(i)
			# ==================================================================
		present.sort(key = sortLast)
		STALA2 = TTT/sqrt(j + 10)
		for i in range(int(sqrt(TRIES)), TRIES):
			for l in range(0, 4):
				if l != 0 or LC0:
					present[i][l] = present[9 - int(floor(sqrt(i)))][l] * (1 + STALA2*(1 - 2*random.random()))
		print(present[0][4])
		#print("==============================================")
		#for j in range(0, 5):
		#	print(present[0][j])
	#var4.set("ER_GEN = " + str(distance_generator()))
	var5.set("ER_OPT = " + str(distance_present2(0)))

# =========================================================================================================================================

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
		if TRYB:
			y = 100*value_present1(0, x)
		else:
			y = 100*value_present2(0, x)
		# MAIN PART
		return 670 - y

def paint():
	global MoveX, MoveY, line_new, cord, X_axis, Y_axis, present
	C.delete(ALL)
	if(len(present[0]) < 2):
		present[0] = [random.random() for i in range(0, 4)]
	#print(optimum[0])
	for w in range(0, 114):
		cord = [funkcja_opt(i) for i in range(w*12, (w + 1)*12 + 2)]
		line_new = C.create_line(cord)
	for i in range(0, len(px)):
		#print(str(px[i]) + "  " + str(py[i]))
		C.create_oval(px[i]*1366, 670 - py[i]*100, px[i]*1366, 670 - py[i]*100, width = 3, fill = 'red')
	#for s in range(0, 4):
	#	for x in range(0, TRIES):
	#		present[0][s] = present[x][s] + present[0][s]
	#	present[0][s] = present[0][s] / TRIES
	#for i in range(0, TRIES):
	#	print("=============  WYNIK  ==================")
	#	print(present[i][0])
	#	print(present[i][1])
	#	print(present[i][2])
	#	print(present[i][3])
	#	print("blad:")
	#	print(present[i][4])
	print("------------------------------------------------------------------")
	print(present[0][0])
	print(present[0][1])
	print(present[0][2])
	print(present[0][3])
	print("opt blad:")
	print(present[0][4])
	X_axis = C.create_line(0, 670, 1366, 670, fill = 'red')
	Y_axis = C.create_line(2, 0, 2, 670, fill = 'red')


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
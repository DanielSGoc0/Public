
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
tabx = []
taby = []

present = []

C = Tkinter.Canvas(root, height=RozmiarY, width=RozmiarX, bg = 'white')
C.pack()

line_old = C.create_line(cord, fill = 'white')
line_new = C.create_line(cord)

X_axis = C.create_line(cord, fill = 'red')
Y_axis = C.create_line(cord, fill = 'red')

#================================================   KONWERSJA  ==========================================================

def convert():
	global tabx, taby
	tabx = [0 for i in range(0, DATA_SIZE)]
	taby = [0 for i in range(0, DATA_SIZE)]
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
				tabx[ind] = x
			else:
				taby[ind] = x
				number = 0
				ind = ind + 1
			number = number + 1
			x = 0
			dot = -1
	tab0 = [i for i in range(0, DATA_SIZE)]
	random.shuffle(tab0)
	px = [0 for i in range(0, DATA_SIZE)]
	py = [0 for i in range(0, DATA_SIZE)]
	for i in range(0, N):
		px[i] = tabx[tab0[i]]
		py[i] = taby[tab0[i]]


#=============================================================================================================================

def value_present(k, x):
	global present
	return present[k][0] * (present[k][1] / (x / present[k][2] + 1) + present[k][3])

def distance_present(k):
	global px, py
	D = 0.0
	for i in range(0, N):
		dy = value_present(k, px[i]) - py[i]
		D += dy*dy
	return D

def da(k):
	global present, px, py, N
	DA = 0.0
	for i in range(0, N):
		DA += (value_present(k, px[i]) - py[i]) * (present[k][1] / (px[i] / present[k][2] + 1) + present[k][3])
	return 2*DA

def db(k):
	global present, px, py, N
	DB = 0.0
	for i in range(0, N):
		DB += (value_present(k, px[i]) - py[i]) * present[k][0] / (px[i] / present[k][2] + 1)
	return 2*DB

def dc(k):
	global present, px, py, N
	DC = 0.0
	for i in range(0, N):
		DC += (value_present(k, px[i]) - py[i]) * present[k][0] * present[k][1] * px[i] / ((px[i] + present[k][2]) * (px[i] + present[k][2]))
	return 2*DC

def dd(k):
	global present, px, py, N
	DD = 0.0
	for i in range(0, N):
		DD += (value_present(k, px[i]) - py[i]) * present[k][0]
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
	for j in range(0, 2**AMT):
		STALA = SSS/sqrt(j + 10)
		for i in range(0, TRIES):
			# ulepszanie
			Da = da(i)
			Db = db(i)
			Dc = dc(i)
			Dd = dd(i)
			if LC0:
				if fabs(present[i][0]) * STALA < fabs(Da):
					if Da > 0:
						present[i][0] = present[i][0] * (1 - 2*random.random() * STALA)
					else:
						present[i][0] = present[i][0] * (1 + 2*random.random() * STALA)
				else:
					present[i][0] -= Da * STALA
			if fabs(present[i][1]) * STALA < fabs(Db):
				if Db/present[i][1] > 0:
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
			present[i][4] = distance_present(i)
			# ==================================================================
		present.sort(key = sortLast)
		STALA2 = TTT/sqrt(j + 10)
		for i in range(int(sqrt(TRIES)), TRIES):
			for l in range(0, 4):
				if l != 0 or LC0:
					present[i][l] = present[int(ceil(sqrt(TRIES))) - int(floor(sqrt(i)))][l] * (1 + STALA2*(1 - 2*random.random()))
	var3.set("ER_OPT = " + str(distance_present(0)))

# =========================================================================================================================================

def Action1():
	global N
	if N > 3:
		N -= 1
		var1.set("N = " + str(N))
		paint()

def Action2():
	global N
	if N < DATA_SIZE:
		N += 1
		var1.set("N = " + str(N))
		paint()

def Action3():
	global AMT
	if AMT > 0:
		AMT -= 1
		var2.set("AMT = " + str(AMT))
		paint()

def Action4():
	global AMT
	if AMT < 15:
		AMT += 1
		var2.set("AMT = " + str(AMT))
		paint()

def Action5():
	convert()
	find_optimum()
	paint()

def Action6():
	global AMT
	for i in range(0, 20):
		convert()
		find_optimum()
		paint()
		tab = [round(-log10(px[i])) for i in range(0, N)]
		tab.sort()
		s = ""
		s = s + str(-round(log10(present[0][2]), 2))
		for i in range(0, N):
			s = s + "   "
			s = s + str(tab[i])
		print(s)

def funkcja_opt(x):
	if x % 2 == 0:
		return x
	else:
		x /= 1366.0
		# MAIN PART
		# argument to obecne x
		y = 100*value_present(0, x)
		# MAIN PART
		return 670 - y

def print_present(k):
	print("------------------------------------------------------------------")
	print(present[k][0])
	print(present[k][1])
	print(present[k][2])
	print(present[k][3])
	print("opt blad:")
	print(present[k][4])

def paint():
	global MoveX, MoveY, line_new, cord, X_axis, Y_axis, present
	C.delete(ALL)
	if(len(present[0]) < 2):	
		present[0] = [random.random() for i in range(0, 4)]
	#print(optimum[0])
	for w in range(0, 114):
		cord = [funkcja_opt(i) for i in range(w*12, (w + 1)*12 + 2)]
		line_new = C.create_line(cord)
	for i in range(0, N):
		C.create_oval(px[i]*1366, 670 - py[i]*100, px[i]*1366, 670 - py[i]*100, width = 3, fill = 'red')
	#print_present(0)
	X_axis = C.create_line(0, 670, 1366, 670, fill = 'red')
	Y_axis = C.create_line(2, 0, 2, 670, fill = 'red')


Button1 = Tkinter.Button(frame, text ="N--", command = Action1)
Button1.pack(side = LEFT)

Button2 = Tkinter.Button(frame, text ="N++", command = Action2)
Button2.pack(side = LEFT)

Button3 = Tkinter.Button(frame, text ="AMT--", command = Action3)
Button3.pack(side = LEFT)

Button3 = Tkinter.Button(frame, text ="AMT++", command = Action4)
Button3.pack(side = LEFT)

Button4 = Tkinter.Button(frame, text ="data", command = Action5)
Button4.pack(side = LEFT)

Button5 = Tkinter.Button(frame, text ="magia...", command = Action6)
Button5.pack(side = LEFT)

Label1 = Label(frame, textvariable=var1, relief=RAISED)
var1.set("N = 13")
Label1.pack(side = LEFT)

Label2 = Label(frame, textvariable=var2, relief=RAISED)
var2.set("AMT = 10")
Label2.pack(side = LEFT)

Label3 = Label(frame, textvariable=var3, relief=RAISED)
var3.set("ER_OPT = ?")
Label3.pack(side = LEFT)

root.mainloop()
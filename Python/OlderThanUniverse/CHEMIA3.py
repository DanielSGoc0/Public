
from Tkinter import *
import Tkinter
from math import *
import random

RozmiarX = 1366
RozmiarY = 670

MoveX = RozmiarX/2
MoveY = RozmiarY/2

DEBUG = FALSE

TRIES = 100
AMT = 0
N = 1
DATA_SIZE = 13
ZERO_MIN = -500 

root = Tk()
frame = Frame(root)
frame.pack()

var1 = StringVar()
var2 = StringVar()
cord = [0, 0, 0, 0]

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

class NUM:
	e = 0
	v = 1

	def __init__(self, exp, val):
		self.e = exp
		self.v = val

class APPRO:
	Ka = [NUM(0, 1) for i in range(0, N + 1)]
	e = [NUM(0, 1) for i in range(0, N + 1)]
	DKa = [NUM(0, 1) for i in range(0, N + 1)]
	De = [NUM(0, 1) for i in range(0, N + 1)]
	ERR = 0
	def __init__(self, X, Y):
		self.Ka = [X[i] for i in range(0, N + 1)]
		self.e = [Y[i] for i in range(0, N + 1)]

def sortNUM(val): 
	return -LOG(val)

def rand_APPRO():
	X = [NUM(random.randint(-40, 20), random.random() + 1) for i in range(0, N + 1)]
	X.sort(key = sortNUM)
	X[0] = NUM(0, 1.0)
	for i in range(1, N + 1):
		X[i] = multiply(X[i], X[i - 1])
	Y = [NUM(random.randint(6, 13), random.random() + 1) for i in range(0, N + 1)]
	return APPRO(X, Y)

#==========================================  ARYTMETYKA  ===========================================

def normalize(C):
	global ZERO_MIN
	while abs(C.v) < 1:
		C.v *= 2
		C.e -= 1
		if C.e < ZERO_MIN:
			C.v = 1
			break
	while abs(C.v) >= 2:
		C.v /= 2
		C.e += 1
	return C

def add(A, B):
	C = NUM(0, 1)
	if A.e > B.e:
		if A.e - B.e > 30:
			C = A
		else:
			C.v = A.v + (B.v / 2**(A.e - B.e))
			C.e = A.e
	else:
		if B.e - A.e > 30:
			C = B
		else:
			C.v = (A.v / 2**(B.e - A.e)) + B.v
			C.e = B.e
	return normalize(C)

def subtract(A, B):
	C = NUM(0, 1)
	if A.e > B.e:
		if A.e - B.e > 30:
			C = A
		else:
			C.v = A.v - B.v / 2**(A.e - B.e)
			C.e = A.e
	else:
		if B.e - A.e > 30:
			C = B
		else:
			C.v = A.v / 2**(B.e - A.e) - B.v
			C.e = B.e
	return normalize(C)

def multiply(A, B):
	C = NUM(A.e + B.e, A.v * B.v)
	return normalize(C)

def divide(A, B):
	C = NUM(A.e - B.e, A.v / B.v)
	return normalize(C)

def ABS(A):
	A.v = abs(A.v)
	return A

def power(A, p):
	if p == 0:
		return NUM(0, 1)
	if p < 0:
		A = divide(NUM(0, 1), A)
		p = -p
	C = A
	for i in range(1, p):
		C = multiply(C, A)
	return C

def printNUM(A):
	print(str(A.v) + " * 2^" + str(A.e))

def strNUM(A):
	return str(A.v) + " * 2^" + str(A.e)

def LOG(A):
	return A.e + log(abs(A.v), 2)

def SIGN(A):
	if(A.v > 0):
		return 1
	return -1

#================================================   KONWERSJA  ==========================================================

def convert():
	global tabx, taby, px, py
	tabx = [NUM(0, 1) for i in range(0, DATA_SIZE)]
	taby = [NUM(0, 1) for i in range(0, DATA_SIZE)]
	plik = open('Chemia/IN3_BM.txt')
	s = plik.read()
	#print s
	number = 1
	ind = 0
	dot = -1
	x = 0.0
	START = False
	minus = True
	for i in range(0, len(s)):
		if 0 <= ord(s[i]) - 48 and ord(s[i]) - 48 <= 9:
			if dot >= 0:
				dot += 1
			if START:
				x *= 10
				x += ord(s[i]) - 48
			else:
				START = True
				x = ord(s[i]) - 48
		elif s[i] == '-':
			minus = True
		elif s[i] == '.' or s[i] == ',':
			dot = 0
		elif START:
			START = False
			x *= 1.0000000000001
			if minus:
				x = -x
			minus = False
			W = normalize(NUM(0, x))
			while dot > 0:
				W.v /= 10
				dot -= 1
				while abs(W.v) < 1:
					W.v *= 2
					W.e -= 1 
			if number == 1:
				tabx[ind] = W
			else:
				taby[ind] = W
				number = 0
				ind += 1
			number += 1
			x = 0
			dot = -1
	#C0 = normalize(NUM(0, 0.0000562))       # CE T
	C0 = normalize(NUM(0, 0.0000113636))   # BM
	px = [NUM(0, 1) for i in range(0, DATA_SIZE)]
	py = [NUM(0, 1) for i in range(0, DATA_SIZE)]
	for i in range(0, DATA_SIZE):
		px[i] = tabx[i]
		py[i] = divide(taby[i], C0)


#=============================================================================================================================


def L_present(k, x):
	global present
	L = NUM(ZERO_MIN, 1.0)
	for i in range(0, N + 1):
		L = multiply(L, x)
		L = add(L, multiply(present[k].e[i], present[k].Ka[i]))
	return L

def M_present(k, x):
	global present
	M = NUM(ZERO_MIN, 1.0)
	for i in range(0, N + 1):
		M = multiply(M, x)
		M = add(M, present[k].Ka[i])
	return M

def distance_present(k):
	global px, py
	D = NUM(ZERO_MIN, 1.0)
	for i in range(0, DATA_SIZE):
		dy = subtract(divide(L_present(k, px[i]), M_present(k, px[i])), py[i])
		D = add(D, multiply(dy, dy))
	return D

def sortLast(val): 
	return LOG(val.ERR)

def FABRICATE():
	global present, px, py
	if len(present) < 2:
		present = [rand_APPRO() for i in range(0, 1)]
	if len(px) < DATA_SIZE or len(py) < DATA_SIZE:
		px = [NUM(0, 1) for i in range(0, DATA_SIZE)]
		py = [NUM(0, 1) for i in range(0, DATA_SIZE)]
	T = [NUM(0, 1.0) for i in range(0, DATA_SIZE)]
	while TRUE:
		present[0] = rand_APPRO()
		if CHECK(0) == TRUE:
			H = NUM(0, 1.0)
			B = TRUE
			for k in range(0, DATA_SIZE):
				px[k] = H
				T[k] = divide(L_present(0, H), M_present(0, H))
				if T[k].e < 6 or T[k].e > 17:
					B = FALSE
					break
				H = divide(H, NUM(0, 10.0))
			if B:
				break
	for k in range(0, DATA_SIZE):
		py[k] = add(T[k], NUM(0, 10.0 - 3 * random.randint(0, 6)))
	present[0].ERR = distance_present(0)
	if DEBUG:
		print("FABRICATED:")
		print_present(0)
		print("END FABRICATED")
		

# ==============================================================  OPTIMUM  ===============================================================================

def CHECK(i):
	global present, N
	if present[i].Ka[1].e >= 24:
		return FALSE
	if N >= 1:
		if present[i].Ka[N].e - present[i].Ka[N - 1].e < -47:
			return FALSE
	for j in range(2, N + 1):
		if divide(multiply(present[i].Ka[j - 1], present[i].Ka[j - 1]), multiply(present[i].Ka[j], present[i].Ka[j - 2])).e < 0:
			return FALSE
	for j in range(0, N + 1):
		if present[i].e[j].e < 3 or present[i].e[j].e > 20:
			return FALSE
	return TRUE

def UPDATE(i):
	global present, TRIES, px, py, ZERO_MIN
	L = [L_present(i, px[t]) for t in range(0, DATA_SIZE)]
	M = [M_present(i, px[t]) for t in range(0, DATA_SIZE)]
	S1Ka = [NUM(ZERO_MIN, 1.0) for t in range(0, DATA_SIZE)]
	S2Ka = [NUM(ZERO_MIN, 1.0) for t in range(0, DATA_SIZE)]
	S1e = [NUM(ZERO_MIN, 1.0) for t in range(0, DATA_SIZE)]
	S2e = [NUM(ZERO_MIN, 1.0) for t in range(0, DATA_SIZE)]
	for k in range(0, DATA_SIZE):
		F = divide(subtract(divide(L[k], M[k]), py[k]), M[k])
		F.e += 1
		for s in range(0, N + 1):
			W = multiply(multiply(F, present[i].Ka[N - s]), subtract(present[i].e[N - s], divide(L[k], M[k])))
			S1Ka[N - s] = add(S1Ka[N - s], W)
			S2Ka[N - s] = add(S2Ka[N - s], ABS(W))
			W = multiply(multiply(F, present[i].e[N - s]), present[i].Ka[N - s])
			S1e[N - s] = add(S1e[N - s], W)
			S2e[N - s] = add(S2e[N - s], ABS(W))
			F = multiply(F, px[k])
	for s in range(0, N + 1):
		present[i].DKa[s] = LOG(add(NUM(0, 1.0), divide(multiply(S1Ka[s], S1Ka[s]), multiply(S2Ka[s], NUM(0, DATA_SIZE))))) * SIGN(S1Ka[s])
		present[i].De[s] = LOG(add(NUM(0, 1.0), divide(multiply(S1e[s], S1e[s]), multiply(S2e[s], NUM(0, DATA_SIZE))))) * SIGN(S1e[s])
		if DEBUG:
			print("s = " + str(s) + "   val DKa: " + str(present[i].DKa[s]) + "   val De: " + str(present[i].De[s]))

def STEP(STALA1, STALA2):
	global present, TRIES, px, py
	for i in range(0, TRIES):
		present[i].ERR = distance_present(i)
	present.sort(key = sortLast)
	if DEBUG:
		print("=========================  ULEPSZANIE  =======================")
	for i in range(0, int(sqrt(TRIES))):
		UPDATE(i)
	for i in range(int(sqrt(TRIES - 1)), TRIES):
		w = int(floor(TRIES/((float) (i + 1)))) - 1
		if DEBUG:
			print("biore " + str(i) + " i ulepszam na " + str(w))
			print("biore " + str(i) + " i ulepszam na " + str(w))
			print("ERR: " + strNUM(present[w].ERR))
			print("LOGLOG: " + str(LOG(add(present[w].ERR, NUM(0, 1.0)))))
		mutation1 = LOG(add(present[w].ERR, NUM(0, 1.0))) * STALA1
		mutation2 = LOG(add(present[w].ERR, NUM(0, 1.0))) * STALA2
		if DEBUG:
			print("mutation1: " + str(mutation1))
			print("mutation2: " + str(mutation2))
		for j in range(0, N + 1):
			if DEBUG:
				print("zmiana j = " + str(j) + ": " + str(-mutation1 * present[w].DKa[j]))
			if j > 0:
				present[i].Ka[j] = multiply(present[w].Ka[j], NUM(0, 2**(-mutation1 * random.random() * present[w].DKa[j])))
			present[i].e[j] = multiply(present[w].e[j], NUM(0, 2**(-mutation2 * random.random() * present[w].De[j])))
	for i in range(0, TRIES):
		if CHECK(i) == FALSE:
			present[i] = rand_APPRO()

	if DEBUG:
		print_present(0)
	#for i in range(int(sqrt(TRIES - 1)), TRIES):
	#	for j in range(0, N):
	#		w = int(ceil(sqrt(TRIES/((float) (i + 1))))) - 1
	#		mutation1 = LOG(add(present[w].ERR, NUM(0, 1))) * random.random()/2
	#		mutation2 = LOG(add(present[w].ERR, NUM(0, 1))) * random.random()/5
	#		present[i].Ka[j] = add(present[w].Ka[j], normalize(NUM(0, mutation1 * present[w].DKa[j])))
	#		present[i].e[j] = add(present[w].e[j], normalize(NUM(0, mutation2 * present[w].De[j])))

def find_optimum():
	global present, TRIES, px, py
	present = [rand_APPRO() for i in range(0, TRIES)]
	if DEBUG:
		print("GENERATED")
		print_present(0)
		print("END GENERATED")
	for j in range(0, 2**AMT):
		#STEP(1.0/(20*log(20 + 2*j)), 1.0/(40*log(60 + 2*j)))
		STEP(1.0/(10*sqrt(20 + 2*j)), 1.0/(40*sqrt(60 + 2*j)))

# =========================================================================================================================================

def Action1():
	global AMT
	if AMT > 0:
		AMT -= 1
		var2.set("AMT = " + str(AMT))
		paint()

def Action2():
	global AMT
	if AMT < 15:
		AMT += 1
		var2.set("AMT = " + str(AMT))
		paint()

def Action3():
	convert()
	find_optimum()
	paint()

def Action4():
	FABRICATE()
	paint()

def Action5():
	find_optimum()
	paint()

def Action6():
	global AMT
	for i in range(0, 2**AMT):
		STEP(1.0/1000, 1.0/2000)
	paint()

def Action7():
	print_present(0)

def funkcja_opt(x):
	if x % 2 == 0:
		return x
	else:
		x /= 1366.0
		# MAIN PART
		# argument to obecne x
		W = normalize(NUM(0, 10**(-x*14.0)))
		S = divide(L_present(0, W), M_present(0, W))
		# MAIN PART
		return 600 - 30.0*LOG(S)

def print_present(k):
	print("------------------------------------------------------------------")
	for i in range(1, N + 1):
		print("Ka " + str(i) + ": " + str(-LOG(divide(present[k].Ka[i], present[k].Ka[i - 1])) * 0.301029996))
	print("--------------")
	for i in range(0, N + 1):
		print("e " + str(i) + ": " + strNUM(present[k].e[i]))
	print("ERR:")
	printNUM(distance_present(k))

def paint():
	global MoveX, MoveY, line_new, cord, X_axis, Y_axis, present, px, py
	C.delete(ALL)
	if(len(present) < TRIES):
		present = [rand_APPRO() for i in range(0, TRIES)]
	for w in range(0, 114):
		cord = [funkcja_opt(i) for i in range(w*12, (w + 1)*12 + 2)]
		line_new = C.create_line(cord)
	for i in range(0, DATA_SIZE):
		C.create_oval(-(LOG(px[i])*1366.0)/47, 600.0 - LOG(py[i])*30.0, -(LOG(px[i])*1366.0)/47, 600.0 - LOG(py[i])*30.0, width = 3, fill = 'red')
	#print_present(0)
	X_axis = C.create_line(0, 670, 1366, 670, fill = 'red')
	Y_axis = C.create_line(2, 0, 2, 670, fill = 'red')


Button1 = Tkinter.Button(frame, text ="AMT--", command = Action1)
Button1.pack(side = LEFT)

Button2 = Tkinter.Button(frame, text ="AMT++", command = Action2)
Button2.pack(side = LEFT)

Button3 = Tkinter.Button(frame, text ="data", command = Action3)
Button3.pack(side = LEFT)

Button4 = Tkinter.Button(frame, text ="FABRICATE", command = Action4)
Button4.pack(side = LEFT)

Button5 = Tkinter.Button(frame, text ="GENERATE", command = Action5)
Button5.pack(side = LEFT)

Button6 = Tkinter.Button(frame, text ="STEPS", command = Action6)
Button6.pack(side = LEFT)

Button7 = Tkinter.Button(frame, text ="PRINT", command = Action7)
Button7.pack(side = LEFT)

Label1 = Label(frame, textvariable=var1, relief=RAISED)
var1.set("N = " + str(N))
Label1.pack(side = LEFT)

Label2 = Label(frame, textvariable=var2, relief=RAISED)
var2.set("AMT = " + str(AMT))
Label2.pack(side = LEFT)

root.mainloop()
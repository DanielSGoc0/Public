
from Tkinter import *
import Tkinter
from math import *
import random

RozmiarX = 1366
RozmiarY = 670

MoveX = RozmiarX/2
MoveY = RozmiarY/2

DEBUG = FALSE

TRIES = 1
AMT = 0
N = 2
DATA_SIZE = 14
ZERO_MIN = -100 

LAMBDA = 615
LAMBDA_UP = 899
LAMBDA_DOWN = 400

root = Tk()
frame = Frame(root)
frame.pack()

var1 = StringVar()
var2 = StringVar()
cord = [0, 0, 0, 0]

TAB = []
px = []

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
	DKa = [NUM(0, 1) for i in range(0, N + 1)]
	e = [NUM(0, 1) for i in range(0, N + 1)]
	ERR = 0
	def __init__(self, X, Y):
		self.Ka = [X[i] for i in range(0, N + 1)]
		self.e = [Y[i] for i in range(0, N + 1)]

def sortNUM(val): 
	return -LOG(val)

def rand_APPRO():
	X = [NUM(random.randint(-40, 0), random.random() + 1) for i in range(0, N + 1)]
	X.sort(key = sortNUM)
	X[0] = NUM(0, 1.0)
	for i in range(1, N + 1):
		X[i] = mul(X[i], X[i - 1])
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

def sub(A, B):
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

def mul(A, B):
	C = NUM(A.e + B.e, A.v * B.v)
	return normalize(C)

def div(A, B):
	C = NUM(A.e - B.e, A.v / B.v)
	return normalize(C)

def ABS(A):
	A.v = abs(A.v)
	return A

def power(A, p):
	if p == 0:
		return NUM(0, 1)
	if p < 0:
		A = div(NUM(0, 1), A)
		p = -p
	C = A
	for i in range(1, p):
		C = mul(C, A)
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
	global TAB, px, DATA_SIZE, WYNIKI
	TAB = [[NUM(0, 1) for i in range(0, DATA_SIZE)] for j in range(0, LAMBDA_UP - LAMBDA_DOWN + 1)]
	px = [NUM(0, 1) for j in range(0, DATA_SIZE)]
	WYNIKI = [rand_APPRO() for j in range(0, 500)]
	#print("===============================   WYGENEROWANE   ======================================")
	plik = open('Chemia/IN4_CET.txt')
	s = plik.read()
	#print s
	ind = -DATA_SIZE
	dot = -1
	x = 0
	START = False
	minus = True
	C0 = normalize(NUM(0, 0.0000562))	 # CE T
	#C0 = normalize(NUM(0, 0.0000113636))   # BM
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
			if W.v < 0:
				W = NUM(ZERO_MIN, 1.0)
			if ind < 0:
				px[DATA_SIZE + ind] = W
			else:
				TAB[ind // DATA_SIZE][ind % DATA_SIZE] = div(W, C0)
			ind += 1
			x = 0
			dot = -1

# ====================================== FABRICATE ====================================================================

def WORKS():
	global present, px
	if CHECK(0):
		B = TRUE
		for i in range(0, DATA_SIZE):
			T = div(add(mul(px[i], present[0].e[0]), mul(present[0].Ka[1], present[0].e[1])), add(px[i], present[0].Ka[1]))
			if T.e < 4 or T.e > 17:
				B = FALSE
				break
		return B
	return FALSE

def FABRICATE():
	global present, TAB
	if len(present) < 2:
		present = [rand_APPRO() for i in range(0, 1)]
	T = [NUM(0, 1.0) for i in range(0, DATA_SIZE)]
	while TRUE:
		present[0] = rand_APPRO()
		if WORKS():
			break
	for i in range(LAMBDA_DOWN, LAMBDA_UP + 1):
		while TRUE:
			present[0].e = [NUM(random.randint(6, 13), random.random() + 1) for j in range(0, N + 1)]
			if WORKS():
				for k in range(0, DATA_SIZE):
					T = div(add(mul(px[k], present[0].e[0]), mul(present[0].Ka[1], present[0].e[1])), add(px[k], present[0].Ka[1]))
					TAB[LAMBDA_UP - i][k] = add(T, NUM(0, 10.0 - 3 * random.randint(0, 6)))
				break

	

#=============================================================================================================================

def UPDATE(k, s, scalar):
	global TAB, px, present

	S0 = NUM(ZERO_MIN, 1.0)
	S1 = NUM(ZERO_MIN, 1.0)
	A = NUM(ZERO_MIN, 1.0)
	B = NUM(ZERO_MIN, 1.0)
	C = NUM(ZERO_MIN, 1.0)
	DS0 = NUM(ZERO_MIN, 1.0)
	DS1 = NUM(ZERO_MIN, 1.0)
	DA = NUM(ZERO_MIN, 1.0)
	DB = NUM(ZERO_MIN, 1.0)
	DC = NUM(ZERO_MIN, 1.0)

	M = [NUM(0, 1.0) for i in range(0, DATA_SIZE)]
	F = [[NUM(0, 1.0) for i in range(0, DATA_SIZE)] for j in range(0, 2)]
	DF = [[NUM(0, 1.0) for i in range(0, DATA_SIZE)] for j in range(0, 2)]

	for i in range(0, DATA_SIZE):
		M[i] = add(px[i], present[k].Ka[1])
		F[0][i] = div(px[i], M[i])
		F[1][i] = div(present[k].Ka[1], M[i])
		DF[0][i] = mul(div(mul(px[i], present[k].Ka[1]), mul(M[i], M[i])), NUM(0, -1.0))
		DF[1][i] = div(mul(px[i], present[k].Ka[1]), mul(M[i], M[i]))

	for i in range(0, DATA_SIZE):
		S0 = add(S0, mul(TAB[s][i], F[0][i]))
		S1 = add(S1, mul(TAB[s][i], F[1][i]))
		A = add(A, mul(F[0][i], F[0][i]))
		B = add(B, mul(F[0][i], F[1][i]))
		C = add(C, mul(F[1][i], F[1][i]))
		DS0 = add(DS0, mul(TAB[s][i], DF[0][i]))
		DS1 = add(DS1, mul(TAB[s][i], DF[1][i]))
		DA = add(DA, mul(mul(DF[0][i], F[0][i]), NUM(1, 1.0)))
		DB = add(DB, add(mul(DF[0][i], F[1][i]), mul(F[0][i], DF[1][i])))
		DC = add(DC, mul(mul(DF[1][i], F[1][i]), NUM(1, 1.0)))
	
	E = [NUM(0, 1.0) for i in range(0, 2)]
	DE = [NUM(0, 1.0) for i in range(0, 2)]
	ME = sub(mul(A, C), mul(B, B))
	DME = sub(add(mul(DA, C), mul(A, DC)), mul(mul(B, DB), NUM(1, 1.0)))

	LE0 = sub(mul(S0, C), mul(S1, B))
	LE1 = sub(mul(A, S1), mul(B, S0))

	E[0] = div(LE0, ME)
	E[1] = div(LE1, ME)

	DE[0] = mul(sub(add(mul(DS0, C), mul(S0, DC)), add(mul(DS1, B), mul(S1, DB))), ME)
	DE[0] = div(sub(DE[0], mul(DME, LE0)), mul(ME, ME))

	DE[1] = mul(sub(add(mul(DA, S1), mul(A, DS1)), add(mul(DB, S0), mul(B, DS0))), ME)
	DE[1] = div(sub(DE[1], mul(DME, LE1)), mul(ME, ME))

	for i in range(0, DATA_SIZE):
		W = add(add(mul(DE[0], F[0][i]), mul(E[0], DF[0][i])), add(mul(DE[1], F[1][i]), mul(E[1], DF[1][i])))
		W = mul(mul(W, NUM(1, 1.0)), sub(add(mul(E[0], F[0][i]), mul(E[1], F[1][i])), TAB[s][i]))
		present[k].DKa[1] = add(present[k].DKa[1], mul(W, scalar))

	present[k].ERR = NUM(ZERO_MIN, 1.0)
	for i in range(0, DATA_SIZE):
		W = sub(add(mul(E[0], F[0][i]), mul(E[1], F[1][i])), TAB[s][i])
		present[k].ERR = add(present[k].ERR, mul(W, W))

	for i in range(0, 2):
		present[k].e[i] = E[i]

def UPDATE_ALL(k):
	global present
	present[k].DKa[1] = NUM(ZERO_MIN, 1.0)
	#UPDATE(k, LAMBDA_UP - LAMBDA, NUM(0, 1.0))
	for i in range(0, LAMBDA_UP - LAMBDA_DOWN + 1):
		UPDATE(k, i, NUM(0, 1/1000.0))

def sortLast(val): 
	return LOG(val.ERR)

# ==============================================================  OPTIMUM  ===============================================================================

def CHECK(i):
	global present, N
	if present[i].Ka[1].e >= 0:
		return FALSE
	if N >= 1:
		if present[i].Ka[N].e - present[i].Ka[N - 1].e < -47:
			return FALSE
	for j in range(2, N + 1):
		if div(mul(present[i].Ka[j - 1], present[i].Ka[j - 1]), mul(present[i].Ka[j], present[i].Ka[j - 2])).e < 0:
			return FALSE
	return TRUE

def STEP(STALA):
	global present, TRIES, px, py
	for i in range(0, TRIES):
		UPDATE_ALL(i)
	present.sort(key = sortLast)
	print("okokokokokok")
	for i in range(int(sqrt(TRIES - 1)), TRIES):
		w = int(floor(TRIES/((float) (i + 1)))) - 1
		mutation = LOG(add(NUM(0, 1.0), div(mul(present[w].ERR, present[w].ERR), mul(present[w].DKa[1], present[w].DKa[1])))) * STALA * random.random()
		if present[w].DKa[1].v > 0:
			mutation *= -1
		print(strNUM(present[w].ERR))
		print(strNUM(present[w].DKa[1]))
		print(str(mutation))
		present[i].Ka[1] = mul(present[w].Ka[1], NUM(0, 2.0**(mutation)))
	for i in range(0, TRIES):
		if CHECK(i) == FALSE:
			present[i] = rand_APPRO()

def find_optimum():
	global present, TRIES, px, py
	present = [rand_APPRO() for i in range(0, TRIES)]
	for j in range(0, 2**AMT):
		#STEP(1.0/(20*log(20 + 2*j)), 1.0/(40*log(60 + 2*j)))
		STEP(1.0/(sqrt(10 + 2*j)))

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
		STEP(1)
	paint()

def Action7():
	print_present(0)

def Action8():
	global LAMBDA, LAMBDA_UP
	if LAMBDA + 10 <= LAMBDA_UP:
		LAMBDA += 10
	paint()

def Action9():
	global LAMBDA, LAMBDA_DOWN
	if LAMBDA - 10 >= LAMBDA_DOWN:
		LAMBDA -= 10
	paint()

def funkcja_opt(x):
	global present
	if x % 2 == 0:
		return x
	else:
		x /= 1366.0
		# MAIN PART
		# argument to obecne x
		W = normalize(NUM(0, 10**(-x*14.0)))
		S = div(add(mul(W, present[0].e[0]), mul(present[0].Ka[1], present[0].e[1])), add(W, present[0].Ka[1]))
		# MAIN PART
		return 600 - 30.0*LOG(S)

def print_present(k):
	print("------------------------------------------------------------------")
	for i in range(1, N + 1):
		print("Ka " + str(i) + ": " + str(-LOG(div(present[k].Ka[i], present[k].Ka[i - 1])) * 0.301029996))
	print("--------------")
	for i in range(0, N + 1):
		print("e " + str(i) + ": " + strNUM(present[k].e[i]))

def paint():
	global MoveX, MoveY, line_new, cord, X_axis, Y_axis, present, px, LAMBDA_UP, LAMBDA
	C.delete(ALL)
	if(len(present) < TRIES):
		present = [rand_APPRO() for i in range(0, TRIES)]
	UPDATE(0, LAMBDA_UP - LAMBDA, NUM(0, 1.0))
	for w in range(0, 114):
		cord = [funkcja_opt(i) for i in range(w*12, (w + 1)*12 + 2)]
		line_new = C.create_line(cord)
	for i in range(0, DATA_SIZE):
		C.create_oval(-(LOG(px[i])*1366.0)/47, 600.0 - LOG(TAB[LAMBDA_UP - LAMBDA][i])*30.0, -(LOG(px[i])*1366.0)/47, 600.0 - LOG(TAB[LAMBDA_UP - LAMBDA][i])*30.0, width = 3, fill = 'red')
	print_present(0)
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

Button8 = Tkinter.Button(frame, text ="LAMBDA+", command = Action8)
Button8.pack(side = LEFT)

Button9 = Tkinter.Button(frame, text ="LAMBDA-", command = Action9)
Button9.pack(side = LEFT)

Label1 = Label(frame, textvariable=var1, relief=RAISED)
var1.set("N = " + str(N))
Label1.pack(side = LEFT)

Label2 = Label(frame, textvariable=var2, relief=RAISED)
var2.set("AMT = " + str(AMT))
Label2.pack(side = LEFT)

root.mainloop()

from Tkinter import *
import Tkinter
from math import *
import random
import copy

RozmiarX = 1366
RozmiarY = 670

MoveX = RozmiarX/2
MoveY = RozmiarY/2

DEBUG = FALSE

TRIES = 1
AMT = 0
N = 2
DATA_SIZE = 4
ZERO_MIN = -100 

LAMBDA = 615
LAMBDA_UP = 800
LAMBDA_DOWN = 500

SINGLE = FALSE

root = Tk()
frame = Frame(root)
frame.pack()

var1 = StringVar()
var2 = StringVar()
cord = [0, 0, 0, 0]

TAB = []
px = []

present = []
T = []
MAT = []
det = []

WYNIKI = []

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
	X[1] = normalize(NUM(0, 0.000724436))
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
	C = NUM(0, 1.0)
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
	C = NUM(0, 1.0)
	if A.e > B.e:
		if A.e - B.e > 30:
			C = A
		else:
			C.v = A.v - B.v / 2**(A.e - B.e)
			C.e = A.e
	else:
		if B.e - A.e > 30:
			C = MINUS(B)
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

def MINUS(A):
	A.v *= -1
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

def SQRT(A):
	if A.e % 2 == 1:
		A.e -= 1
		A.v *= 2
	A.e /= 2
	A.v = sqrt(A.v)
	return normalize(A)

def VAL(A):
	return 2**A.e * A.v

#================================================   KONWERSJA  ==========================================================


def convert():
	global TAB, px, DATA_SIZE, WYNIKI, LAMBDA_DOWN, LAMBDA_UP, LAMBDA
	TAB = [[NUM(0, 1) for i in range(0, DATA_SIZE)] for j in range(0, LAMBDA_UP - LAMBDA_DOWN + 1)]
	px = [NUM(0, 1) for j in range(0, DATA_SIZE)]
	WYNIKI = [rand_APPRO() for j in range(0, LAMBDA_UP - LAMBDA_DOWN + 1)]
	#print("===============================   WYGENEROWANE   ======================================")
	plik = open('Chemia/anime')
	s = plik.read()
	#print s
	ind = -DATA_SIZE
	dot = -1
	x = 0
	START = False
	minus = False
	#C0 = normalize(NUM(0, 0.0000562))	 # CE T
	C0 = normalize(NUM(0, 0.0000113636))   # BM
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

def VALUE(k, x):
	global present
	L = NUM(ZERO_MIN, 1.0)
	M = NUM(ZERO_MIN, 1.0)
	for i in range(0, N + 1):
		L = mul(L, x)
		L = add(L, mul(present[k].e[i], present[k].Ka[i]))
		M = mul(M, x)
		M = add(M, present[k].Ka[i])
	return div(L, M)

def WORKS():
	global present, px
	if CHECK(0):
		B = TRUE
		for i in range(0, DATA_SIZE):
			T = VALUE(0, px[i])
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
					T = VALUE(0, px[k])
					TAB[LAMBDA_UP - i][k] = add(T, NUM(0, 10.0 - 3 * random.randint(0, 6)))
				break

#=============================================================================================================================

def MATRIX():
	global T, MAT
	RES = NUM(0, 1.0)
	MAT = [[NUM(ZERO_MIN, 1.0) for j in range(0, 2*N + 2)] for i in range(0, N + 1)]
	for i in range(0, N + 1):
		for j in range(0, N + 1):
			MAT[i][j] = T[i][j][0]
	for i in range(0, N + 1):
		MAT[i][i + N + 1] = NUM(0, 1.0)

	for x in range(0, N + 1):
		#for j in range(x + 1, N + 1):
		#	if(LOG(ABS(M[x][j])) > LOG(ABS(M[x][x]))):
		#		SWAP(x, j)
		for j in range(x + 1, 2*N + 2):
			MAT[x][j] = div(MAT[x][j], MAT[x][x])
		RES = mul(RES, MAT[x][x])	
		MAT[x][x] = NUM(0, 1.0)
		for i in range(x + 1, N + 1):
			for j in range(x + 1, 2*N + 2):
				MAT[i][j] = sub(MAT[i][j], mul(MAT[x][j], MAT[i][x]))
			#M[i][x] = NUM(ZERO_MIN, 1.0)
	for x in range(0, N):
		for i in range(0, N - x):
			for j in range(N + 1, 2*N + 2):
				MAT[i][j] = sub(MAT[i][j], mul(MAT[N - x][j], MAT[i][N - x]))
			#M[i][N - x] = NUM(ZERO_MIN, 1.0)
	return RES
	

def UPDATE(k, s, STOP):
	global TAB, px, present, T, dp, POWER, N, MAT, DATA_SIZE, det

	T = [[[NUM(ZERO_MIN, 1.0) for t in range(0, N + 1)] for j in range(0, N + 1)] for i in range(0, N + 1)]
	T_Ka = [[NUM(ZERO_MIN, 1.0) for t in range(0, N + 1)] for j in range(0, N + 1)]
	inv = [[[NUM(ZERO_MIN, 1.0) for t in range(0, N + 1)] for j in range(0, N + 1)] for i in range(0, N + 2)]
	
	M = [NUM(ZERO_MIN, 1.0) for i in range(0, DATA_SIZE)]
	L = [[NUM(ZERO_MIN, 1.0) for i in range(0, N + 1)] for j in range(0, DATA_SIZE)]
	F = [[NUM(ZERO_MIN, 1.0) for i in range(0, N + 1)] for j in range(0, DATA_SIZE)]
	DF = [[[NUM(ZERO_MIN, 1.0) for x in range(0, N + 1)] for i in range(0, N + 1)] for j in range(0, DATA_SIZE)]
	det = [[NUM(ZERO_MIN, 1.0) for i in range(0, N + 1)] for j in range(0, N + 2)]

	for i in range(0, DATA_SIZE):
		Q = NUM(0, 1.0)
		for j in range(0, N + 1):
			L[i][N - j] = mul(Q, present[k].Ka[N - j])
			M[i] = add(M[i], L[i][N - j])
			Q = mul(Q, px[i])

	for i in range(0, DATA_SIZE):
		for j in range(0, N + 1):
			F[i][j] = div(L[i][j], M[i])
			for x in range(0, N + 1):
				if x == j:
					DF[i][j][x] = div(mul(L[i][j], sub(M[i], L[i][j])), mul(M[i], M[i]))
				else:
					DF[i][j][x] = mul(div(mul(L[i][x], L[i][j]), mul(M[i], M[i])), NUM(0, -1.0))
	
	for i in range(0, DATA_SIZE):
		for y in range(0, N + 1):
			T_Ka[y][0] = add(T_Ka[y][0], mul(TAB[s][i], F[i][y]))
			for x in range(0, N + 1):
				T[y][x][0] = add(T[y][x][0], mul(F[i][y], F[i][x]))
				for t in range(1, N + 1):
					T[y][x][t] = add(T[y][x][t], add(mul(DF[i][y][t], F[i][x]), mul(F[i][y], DF[i][x][t])))
			for t in range(1, N + 1):
				T_Ka[y][t] = add(T_Ka[y][t], mul(TAB[s][i], DF[i][y][t]))

	det[N + 1][0] = MATRIX()

	for y in range(0, N + 1):
		for x in range(0, N + 1):
			inv[N + 1][y][x] = MAT[y][N + 1 + x]

	for i in range(0, N + 1):
		for j in range(0, N + 1):
			R = copy.deepcopy(T_Ka[j][i])
			T_Ka[j][i] = copy.deepcopy(T[j][i][0])
			T[j][i][0] = R
		det[i][0] = MATRIX()
		for y in range(0, N + 1):
			for x in range(0, N + 1):
				inv[i][y][x] = MAT[y][N + 1 + x]
		for j in range(0, N + 1):
			R = copy.deepcopy(T_Ka[j][i])
			T_Ka[j][i] = copy.deepcopy(T[j][i][0])
			T[j][i][0] = R

	for i in range(0, N + 1):
		W = NUM(ZERO_MIN, 1.0)
		for j in range(0, N + 1):
			W = add(W, mul(inv[N + 1][i][j], T_Ka[j][0]))
		present[k].e[i] = W

	if STOP:
		return

	# =========================================================

	De = [[NUM(ZERO_MIN, 1.0) for i in range(0, N + 1)] for j in range(0, N + 1)]

	for i in range(0, N + 2):
		for j in range(1, N + 1):
			for y in range(0, N + 1):
				for x in range(0, N + 1):
					if x == i:
						det[i][j] = add(det[i][j], mul(inv[i][y][x], T_Ka[y][j]))
					else:
						det[i][j] = add(det[i][j], mul(inv[i][y][x], T[y][x][j]))
			det[i][j] = div(det[i][j], det[i][0])


	for i in range(0, N + 1):
		for j in range(0, N + 1):
			De[i][j] = div(sub(mul(det[i][j], det[N + 1][0]), mul(det[i][0], det[N + 1][j])), mul(det[N + 1][0], det[N + 1][0]))

	for i in range(0, DATA_SIZE):
		W = NUM(ZERO_MIN, 1.0)
		for j in range(0, N + 1):
			W = add(W, mul(present[k].e[j], F[i][j]))
		W = sub(W, TAB[s][i])
		present[k].ERR = add(present[k].ERR, mul(W, W))
		for j in range(0, N + 1):
			V = NUM(ZERO_MIN, 1.0)
			for x in range(0, N + 1):
				V = add(V, add(mul(De[x][j], F[i][x]), mul(present[k].e[x], DF[i][x][j])))
			V.e += 1
			present[k].DKa[j] = add(present[k].DKa[j], mul(V, W))



def UPDATE_ALL(k):
	global present, LAMBDA, LAMBDA_UP
	for i in range(0, N + 1):
		present[k].DKa[i] = NUM(ZERO_MIN, 1.0)
	present[k].ERR = NUM(ZERO_MIN, 1.0)
	#UPDATE(k, LAMBDA_UP - LAMBDA, NUM(0, 1.0))
	if SINGLE:
		UPDATE(k, LAMBDA_UP - LAMBDA, FALSE)
	else:
		for i in range(0, LAMBDA_UP - LAMBDA_DOWN + 1):
			UPDATE(k, i, FALSE)

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
	for i in range(int(sqrt(TRIES - 1)), TRIES):
		w = int(floor(TRIES/((float) (i + 1)))) - 1
		LL = NUM(ZERO_MIN, 1.0)
		for j in range(0, N + 1):
			LL = add(LL, mul(present[w].DKa[j], present[w].DKa[j]))
		#mutation = - LOG(add(NUM(0, 1.0), div(mul(present[w].ERR, present[w].ERR), LL))) * STALA * random.random()
		mutation = -STALA
		LL = SQRT(LL)
		for j in range(1, N + 1):
			present[i].Ka[j] = mul(present[w].Ka[j], NUM(0, 2.0**(mutation * VAL(div(present[w].DKa[j], LL)))))
	for i in range(0, TRIES):
		if CHECK(i) == FALSE:
			present[i] = rand_APPRO()

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
	global present
	convert()
	present = [rand_APPRO() for i in range(0, TRIES)]
	paint()

def Action4():
	FABRICATE()
	paint()

def Action5():
	global present
	present = [rand_APPRO() for i in range(0, TRIES)]
	paint()

def Action6():
	global AMT
	#STEP(0.5)
	for i in range(0, 2**AMT):
		STEP(2.0/(i + 1))
	paint()

def Action7():
	global WYNIKI
	WYNIKI = [rand_APPRO() for i in range(0, LAMBDA_UP - LAMBDA_DOWN + 1)]
	for i in range(0, LAMBDA_UP - LAMBDA_DOWN + 1):
		UPDATE(0, i, TRUE)
		WYNIKI[i] = copy.deepcopy(present[0])
	print("============================== WYNIKI ===========================")
	for j in range(0, LAMBDA_UP - LAMBDA_DOWN + 1):
		for i in range(0, N):
			WYNIKI[j].Ka[N - i] = div(WYNIKI[j].Ka[N - i], WYNIKI[j].Ka[N - 1 - i])
	for i in range(1, N + 1):
		print("=======  Ka" + str(i) + "  =========")
		for j in range(0, LAMBDA_UP - LAMBDA_DOWN + 1):
			print(str(-LOG(WYNIKI[j].Ka[i]) * 0.301029996))
	for i in range(0, N + 1):
		print("=======  e" + str(i) + "  =========")
		for j in range(0, LAMBDA_UP - LAMBDA_DOWN + 1):
			print(str(VAL(WYNIKI[j].e[i])))
		

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
		S = VALUE(0, W)
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
	UPDATE(0, LAMBDA_UP - LAMBDA, TRUE)
	for w in range(0, 114):
		cord = [funkcja_opt(i) for i in range(w*12, (w + 1)*12 + 2)]
		line_new = C.create_line(cord)
	for i in range(0, DATA_SIZE):
		C.create_oval(-(LOG(px[i])*1366.0)/47, 600.0 - LOG(TAB[LAMBDA_UP - LAMBDA][i])*30.0, -(LOG(px[i])*1366.0)/47, 600.0 - LOG(TAB[LAMBDA_UP - LAMBDA][i])*30.0, width = 5, fill = 'red')
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
from sympy import *
init_printing()

# The idea behind this module is the following:
# Simply, we represent all series as finite lists of Sympy terms.
# We implement the basic series operations by hand.


# ====================  HELPERS AND SERIES  =========================

def print_var(st, s):
	pprint("\n" + st + " = \n")
	pprint(s)

def bits(n):
	res = []
	while n:
		b = n & (~n+1)
		res.append(b)
		n ^= b
	return res

def simpler(terms):
	return [simplify(term) for term in terms]

# returns the series of exponential function
def exp_series(order):
	result = [Integer(0) for k in range(order + 1)]
	result[0] = Integer(1)

	denom = 1
	for k in range(1, order + 1):
		denom *= k
		result[k] = Rational(1, denom)
	return result

# returns the series of exp(c * x**2)
def exp2_series(c, order):
	result = [Integer(0) for k in range(order + 1)]
	result[0] = Integer(1)

	denom = 1
	for k in range(2, order + 1, 2):
		denom *= (k//2)
		result[k] = Rational(1, denom) * c**(k//2)
	return result

def one_series(order):
	result = [Integer(1)]
	for k in range(1, order + 1):
		result.append(Integer(0))
	return result

def zero_series(order):
	return [Integer(0) for k in range(order + 1)]


# ====================  SIMPLIFICATION  =========================

def power(terms, pow, order_of_term, order_of_result):
	result = [Integer(0) for k in range(order_of_result + 1)]
	
	if pow == 0:
		result[0] = Integer(1)
		return result
	elif pow * order_of_term > order_of_result:
		return result

	indices = [order_of_term for i in range(pow)]
	sum = order_of_term * pow

	factorial = 1
	for i in range(1, pow + 1):
		factorial *= i

	while True:
		factor = factorial
		last_index_equal = 0
		for i in range(1, pow):
			if indices[i] == indices[last_index_equal]:
				factor //= (i - last_index_equal + 1)
			elif indices[i] > indices[last_index_equal]:
				last_index_equal = i
			else:
				factor = 0
				break

		if factor != 0:
			new_summand = terms[indices[0]]
			for i in range(1, pow):
				new_summand *= terms[indices[i]]
			result[sum] += Integer(factor) * new_summand

		# calculate next indices sequence
		stop = True
		for i in range(pow):
			indices[i] += 1
			sum += 1
			if sum > order_of_result:
				sum -= indices[i] - order_of_term
				indices[i] = order_of_term
			else:
				stop = False
				break

		if stop:
			break
	
	return result


def power_series(terms1, terms2, order_of_term, order_of_result):
	result = [Integer(0) for k in range(order_of_result + 1)]
	expand_to = order_of_result // order_of_term
	for i in range(0, expand_to + 1):
		if terms1[i] == Integer(0):
			continue

		result_partial = power(terms2, i, order_of_term, order_of_result)
		for j in range(i*order_of_term, order_of_result + 1):
			result[j] += result_partial[j] * terms1[i]
	return result


# calculate e^(terms)
# we assume that order_of_term >= 1
def exponentiate(terms, order_of_term, order_of_result):
	return power_series(exp_series(order_of_result), terms, order_of_term, order_of_result)


# This function might be able to be sped up eventually:
# - Firstly, by assuming that xn and yn are just 1, and then by removing the last row and column
# - Secondly, by precalculating the exponentials, so as to not calculate them multiple times.
# - (and perhaps by applying more row/column operations, e.g. by iterating determinant formulas)
# Calculates the determinant of Aij = exp(xi * yj)
# We assume that termsx[i] has order order_of_termsx[i], and that termsx[i] - termsx[j] has
# order alpha[i][j]. Similarly for termsy.
def determinant(termsx, order_of_termsx, termsy, order_of_termsy, alpha, beta, order_of_result):
	N = len(termsx)
	sum = 0
	for i in range(N):
		for j in range(i + 1, N):
			sum += alpha[i][j]
			sum += beta[i][j]
	sum_alpha = [0 for i in range(N)]
	sum_beta = [0 for j in range(N)]
	for i in range(N):
		for j in range(N):
			if i != j:
				sum_alpha[i] += alpha[i][j]
				sum_beta[j] += beta[j][i]
	
	orders_of_exponentials = [[order_of_result - sum + sum_alpha[i] + sum_beta[j] for j in range(N)] for i in range(N)]

	exponentials = [[[] for j in range(N)] for i in range(N)]
	for i in range(N):
		for j in range(N):
			# print("i", i, "j", j, order_of_termsx[i], order_of_termsx[j], orders_of_exponentials[i][j])
			exponentials[i][j] = multiply(termsx[i], order_of_termsx[i], termsy[j], order_of_termsy[j], orders_of_exponentials[i][j])
			exponentials[i][j] = exponentiate(exponentials[i][j], order_of_termsx[i] + order_of_termsy[j], orders_of_exponentials[i][j])
			for k in range(orders_of_exponentials[i][j], order_of_result):
				exponentials[i][j].append(Integer(0))
	
	DP = [zero_series(order_of_result) for s in range(2**N)]
	DP_order = [0 for s in range(2**N)]

	for s in range(1, 2**N):
		logs2 = [b.bit_length() - 1 for b in bits(s)]
		sub_DP = [s - b for b in bits(s)]
		if len(sub_DP) == 1:
			DP_order[s] = 0
			DP[s] = exponentials[-1][logs2[0]]
			continue
		
		DP_order[s] = DP_order[sub_DP[0]]
		for i in range(1, len(sub_DP)):
			DP_order[s] += alpha[-len(sub_DP)][-i]
			DP_order[s] += beta[logs2[0]][logs2[i]]

		for i in range(len(sub_DP)):
			new_term = multiply(exponentials[-len(sub_DP)][logs2[i]], 0, DP[sub_DP[i]], DP_order[sub_DP[i]], order_of_result)
			new_term = scalar_multiply(new_term, Integer((-1)**i))
			DP[s] = add(DP[s], new_term, order_of_result)

	return DP[2**N - 1]


def det(terms, order_of_terms, order_of_result):
	N = len(terms)
	alpha = [[min(order_of_terms[i], order_of_terms[j]) for j in range(N)] for i in range(N)]
	beta = [[min(order_of_terms[i], order_of_terms[j]) for j in range(N)] for i in range(N)]
	return determinant(terms, order_of_terms, terms, order_of_terms, alpha, beta, order_of_result)


# Divides terms1 by terms2.
def divide(terms1, order1, terms2, order2, order_of_result):
	if order1 < order2:
		raise Exception("we divided and expression by another expression of higher order.")
	
	result = [Integer(0) for k in range(order_of_result + 1)]

	for k in range(order1 - order2, order_of_result + 1):
		result[k] = terms1[k + order2]
		for i in range(k + order2 - order1):
			result[k] -= result[i + order1 - order2] * terms2[2*order2 - order1 + k - i]
		result[k] /= terms2[order2]
		result[k] = cancel(result[k])

	return result


# calculate e^(c * (terms)**2)
# we assume that order_of_term >= 1
def exponentiate_squared(terms, c, order_of_term, order_of_result):
	return power_series(exp2_series(c, order_of_result), terms, order_of_term, order_of_result)


def multiply(terms1, order_of_term1, terms2, order_of_term2, order_of_result):
	result = [Integer(0) for k in range(order_of_result + 1)]

	for s in range(order_of_term1 + order_of_term2, order_of_result + 1):
		for i in range(order_of_term1, s + 1 - order_of_term2):
			result[s] += terms1[i] * terms2[s - i]
	return result


def scalar_multiply(terms, c):
	return [c * term for term in terms]


def add(terms1, terms2, order_of_result):
	return [terms1[k] + terms2[k] for k in range(order_of_result + 1)]


def subtract(terms1, terms2, order_of_result):
	return [terms1[k] - terms2[k] for k in range(order_of_result + 1)]


# we assume that terms[order_of_terms] = t**2
def sqrt_terms(terms, order_of_terms, t, order_of_result):
	if order_of_terms % 2 == 1:
		raise Exception("sqrt of odd power!")
	k = order_of_terms//2
	print_var("APPLYING sqrt_terms. input diff", simplify(terms[order_of_terms] - t**2))

	result = [Integer(0) for k in range(order_of_result + 1)]
	result[k] = t
	for n in range(1, order_of_result + 1 - k):
		result[n + k] = terms[n + 2*k]
		for i in range(1, n):
			result[n + k] -= result[k + i] * result[k + n - i]
		result[n + k] /= (2 * t)
		result[n + k] = cancel(result[n + k])
	return result


# =============================================


a11, a12, a13, a14, a15, a16, a17, a18, a19, a110 = symbols('a11 a12 a13 a14 a15 a16 a17 a18 a19 a110')
a22, a23, a24, a25, a26, a27, a28, a29, a210 = symbols('a22 a23 a24 a25 a26 a27 a28 a29 a210')
a34, a35, a36, a37, a38, a39, a310 = symbols('a34 a35 a36 a37 a38 a39 a310')
z0, z1, z2 = symbols('z0 z1 z2')

a11 = z0
a22 = z0 * z1
a23 = - z0**3 / 2 + z0 * z1**2
a34 = z0 * Rational(1, 4) * (-2 * sqrt(Integer(2)) * a12 * z2 - 7 * z0**2 * z1 + 2 * sqrt(Integer(2)) * z0 * z1 * z2 + 4 * z1**3)
a35 = z0 * Rational(1, 8) * (-8 * sqrt(Integer(2)) * a12 * z1 * z2 - 4 * sqrt(Integer(2)) * a13 * z2 + 5 * z0**4 - 2 * sqrt(Integer(2)) * z0**3 * z2 - 32 * z0**2 * z1**2 + 12 * sqrt(Integer(2)) * z0 * z1**2 * z2 + 8 * z1**4) 
a36 = z0 * Rational(1, 96) * (-24 * a12**2 * z1 + 76 * sqrt(Integer(2)) * a12 * z0**2 * z2 + 48 * a12 * z0 * z1**2 - 48 * a12 * z0 * z2**2 - 144 * sqrt(2) * a12 * z1**2 * z2 - 96 * sqrt(Integer(2)) * a13 * z1 * z2 - 48 * sqrt(Integer(2)) * a14 * z2 + 365 * z0**4 * z1 - 208 * sqrt(2) * z0**3 * z1 * z2 - 744 * z0**2 * z1**3 + 48 * z0**2 * z1 * z2**2 + 288 * sqrt(2) * z0 * z1**3 * z2 + 96 * z1**5)
a37 = z0 * Rational(1, 48) * (-4 * sqrt(Integer(2)) * a12**2 * z0 * z2 - 24 * a12**2 * z1**2 + 24 * a12**2 * z2**2 - 24 * a12 * a13 * z1 - 12 * a12 * z0**3 * z1)
a37 = a37 +	z0 * Rational(1, 48) * (184 * sqrt(Integer(2)) * a12 * z0**2 * z1 * z2 + 72 * a12 * z0 * z1**3 - 144 * a12 * z0 * z1 * z2**2 - 96 * sqrt(Integer(2)) * a12 * z1**3 * z2 + 38 * sqrt(Integer(2)) * a13 * z0**2 * z2)
a37 = a37 +	z0 * Rational(1, 48) * (24 * a13 * z0 * z1**2 - 24 * a13 * z0 * z2**2 - 72 * sqrt(Integer(2)) * a13 * z1**2 * z2 - 48 * sqrt(Integer(2)) * a14 * z1 * z2 - 24 * sqrt(Integer(2)) * a15 * z2 - 49 * z0**6)
a37 = a37 +	z0 * Rational(1, 48) * (34 * sqrt(Integer(2)) * z0**5 * z2 + 666 * z0**4 * z1**2 - 12 * z0**4 * z2**2 - 434 * sqrt(Integer(2)) * z0**3 * z1**2 * z2 - 648 * z0**2 * z1**4 + 144 * z0**2 * z1**2 * z2**2 + 240 * sqrt(Integer(2)) * z0 * z1**4 * z2 + 48 * z1**6)

x0 = [Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)]
x1 = [Integer(0), a11, a12, a13, a14, a15, a16, Integer(0), Integer(0), Integer(0), Integer(0)]
x2 = [Integer(0), a11, a22, a23, a24, Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)]
x3 = [Integer(0), a11, a22, a23, a34, a35, a36, a37, a38, Integer(0), Integer(0)]


det1 = simpler(det([subtract(x0, x3, 9), subtract(x1, x3, 9)], [1, 2], 10))
# print_var("det1", det1)
det2 = simpler(det([subtract(x0, x3, 9), subtract(x1, x3, 10), subtract(x2, x3, 10)], [1, 2, 4], 14))
# print_var("det2", det2)
det3 = simpler(det([subtract(x0, x3, 6), subtract(x1, x3, 8), subtract(x2, x3, 10), zero_series(100)], [1, 2, 4, 10000], 22))
# print_var("det3", det3)

# sum_k Cik * Cjk = exp(-|xi - xj|^2 / 2)
# Dik = Cik * exp(|xi - xN|^2 / 2)
# DNk = CNk
# sum_k Dik * Djk = exp(<xi - xN, xj - xN>)
# Aij = exp(<xi - xN, xj - xN>)

# D00 = exponentiate_squared(subtract(x0, x3, 9), Rational(1, 2), 1, 9)
D00inv = exponentiate_squared(subtract(x0, x3, 10), Rational(-1, 2), 1, 10)

D11 = sqrt_terms(det1, 2, x1[1] - x0[1], 9)
D11 = multiply(D11, 0, D00inv, 0, 9)
D11 = simpler(D11)
print_var("D11", D11)

D22 = divide(det2, 8, det1, 2, 12)
D22 = sqrt_terms(D22, 6, sqrt(Rational(1, 2)) * (x0[1] - x2[1]) * (x1[2] - x2[2]), 9)
D22 = simpler(D22)
print_var("D22", D22)

# D33 = divide(det3, 22, multiply(det2, 8, power(multiply(multiply(subtract(x0, x3, 1), 1, subtract(x1, x3, 2), 2, 3), 3, subtract(x2, x3, 4), 4, 7), 2, 7, 14), 14, 22), 22, 0)
# D33 = sqrt_terms(D33, 0)
# D33 = multiply(D33, 0, multiply(multiply(subtract(x0, x3, 1), 1, subtract(x1, x3, 2), 2, 3), 3, subtract(x2, x3, 4), 4, 7), 7, 7)
# D33 = simpler(D33)
# print_var("D33", D33)

D33squared = divide(det3, 22, det2, 8, 14)
print_var("D33squared", D33squared)


D10 = exponentiate(multiply(subtract(x0, x3, 8), 1, subtract(x1, x3, 9), 2, 10), 3, 10)
D10 = multiply(D10, 0, D00inv, 0, 10)
D10 = simpler(D10)
print_var("D10", D10)

D20 = exponentiate(multiply(subtract(x0, x3, 6), 1, subtract(x2, x3, 9), 4, 10), 5, 10)
D20 = multiply(D20, 0, D00inv, 0, 10)
D20 = simpler(D20)
print_var("D20", D20)

D21 = exponentiate(multiply(subtract(x1, x3, 6), 2, subtract(x2, x3, 8), 4, 10), 6, 10)
D21 = subtract(D21, multiply(D10, 0, D20, 0, 10), 10)
D21 = divide(D21, 2, D11, 1, 9)
D21 = simpler(D21)
print_var("D21", D21)

D30 = D00inv
D30 = simpler(D30)
print_var("D30", D30)

D31 = subtract(one_series(10), multiply(D10, 0, D30, 0, 10), 10)
D31 = divide(D31, 2, D11, 1, 9)
D31 = simpler(D31)
print_var("D31", D31)

D32 = subtract(one_series(10), add(multiply(D20, 0, D30, 0, 10), multiply(D21, 1, D31, 1, 10), 10), 10)
D32 = divide(D32, 6, D22, 3, 7)
D32 = simpler(D32)
print_var("D32", D32)


dmean = add(add(scalar_multiply(D30, z0), scalar_multiply(D31, z1), 7), scalar_multiply(D32, z2), 7)
dmean.insert(0, Integer(0))
dmean = simpler(subtract(x3, dmean, 8))
print_var("dmean", dmean)

D33squared.insert(0, Integer(0))
D33squared.insert(0, Integer(0))
error = add(power(dmean, 2, 8, 16), D33squared, 16)
error = simpler(error)
print_var("error", error) 

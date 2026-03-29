from sympy import *
init_printing()

# The idea behind this module is the following:
# Simply, we represent all series as finite lists of Sympy terms.
# We implement the basic series operations by hand.


# ====================  HELPERS AND SERIES  =========================

def print_var(st, s):
	pprint("\n" + st + " = \n")
	pprint(s)

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

# Calculates the determinant of Aij = exp(xi * yj)
# For efficiency of calculation we assume that order_of_termx and order_of_termy are increasing
def determinant(termsx, order_of_termx, termsy, order_of_termy, order_of_result):
	N = len(termsx)
	if N == 1:
		exponential = multiply2(termsx[0], order_of_termx[0], termsy[k], order_of_termy[0], order_of_result)
		exponential = exponentiate(exponential, order_of_termx[0] + order_of_termy[0], order_of_result)
		return exponential
	result = [Integer(0) for k in range(order_of_result + 1)]

	for k in range(N):
		termsx2 = termsx[1:]
		order_of_termx2 = order_of_termx[1:]
		termsy2 = termsy[:k] + termsy[(k+1):]
		order_of_termy2 = order_of_termy[:k] + order_of_termy[(k+1):]

		order = 0
		for i in range(0, N - 1):
			for j in range(i + 1, N - 1):
				order += order_of_termx2[i]
				order += order_of_termy2[i]

		minor_det = determinant(termsx2, order_of_termx2, termsy2, order_of_termy, order_of_result)
		exponential = multiply2(termsx[0], order_of_termx[0], termsy[k], order_of_termy[0], order_of_result)
		exponential = exponentiate(exponential, order_of_termx[0] + order_of_termy[0], order_of_result)
		minor_det = multiply2(minor_det, order, exponential, 0, order_of_result)
		minor_det = scalar_multiply(minor_det, Integer((-1)**k))
		result = add(result, minor_det, order_of_result)

	result = simpler(result)
	return result



# calculate e^(c * (terms)**2)
# we assume that order_of_term >= 1
def exponentiate_squared(terms, c, order_of_term, order_of_result):
	return power_series(exp2_series(c, order_of_result), terms, order_of_term, order_of_result)


def multiply(terms1, terms2, order_of_result):
	result = [Integer(0) for k in range(order_of_result + 1)]

	for s in range(order_of_result + 1):
		for i in range(s + 1):
			result[s] += terms1[i] * terms2[s - i]
	return result


def multiply2(terms1, order_of_term1, terms2, order_of_term2, order_of_result):
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


# we assume that terms[0] = 1
def inverse(terms, order_of_result):
	result = [Integer(0) for k in range(order_of_result + 1)]
	print_var("APPLYING INVERSE", terms[0])

	result[0] = Integer(1)
	for s in range(order_of_result + 1):
		for i in range(s):
			result[s] -= terms[s - i] * result[i]
	return result


# we assume that terms[0] = 1
def sqrt(terms, order_of_result):
	result = [Integer(0) for k in range(order_of_result + 1)]
	print_var("APPLYING SQRT", terms[0])

	result[0] = Integer(1)
	result[1] = Rational(1, 2) * terms[1]
	for s in range(2, order_of_result + 1):
		result[s] = terms[s]
		for i in range(1, s):
			result[s] -= result[s - i] * result[i]
		result[s] *= Rational(1, 2)
	return result


# =============================================

a11, a12, a13, a14, a15, a16, a17, a18, a19, a110 = symbols('a11 a12 a13 a14 a15 a16 a17 a18 a19 a110')
a22, a23, a24, a25, a26, a27, a28, a29, a210 = symbols('a22 a23 a24 a25 a26 a27 a28 a29 a210')
z0, z1 = symbols('z0 z1')


a11 = z0
a12 = Integer(0)
# a13 = Rational(-1, 2) * z0**3 + 3 * z0
# a14 = Integer(0)
a22 = z0 * z1
a23 = - z0**3 / 2 + z0 * z1**2
a24 = Rational(-7, 4) * z0**3 * z1 + z0 * z1**3
# a25 = Rational(5, 8) * z0**5 - 4 * z0**3 * z1**2 + z0 * z1**4
# a26 = Rational(-1, 2) * z0**3 * z1 + Rational(365, 96) * z0**5 * z1 + Rational(-31, 4) * z0**3 * z1**3 + z0 * z1**5
# d0, d1, d2 = symbols('d0 d1 d2')

x0 = [Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)]
x1 = [Integer(0), a11, a12, a13, a14, a15, a16, a17, a18, Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)]
x2 = [Integer(0), a11, a22, a23, a24, a25, a26, a27, a28, Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)]



D00 = simpler(exponentiate_squared(subtract(x0, x2, 14), Rational(1, 2), 1, 14))
D00inv = simpler(inverse(D00, 14))
print_var("D00inv", D00inv)

D10 = simpler(exponentiate(scalar_multiply(multiply(subtract(x0, x2, 5), subtract(x1, x2, 5), 5), Integer(-1)), 3, 5))
D10 = simpler(multiply(D10, D00inv, 5))
print_var("D10", D10)

D11 = exp2_series(Integer(-1), 14)
D11 = subtract(one_series(14), D11, 14)
D11.pop(0)
D11.pop(0)
D11 = simpler(sqrt(D11, 12))
D11 = power_series(D11, subtract(x1, x0, 12), 1, 12)
D11 = simpler(multiply(D11, exponentiate_squared(subtract(x1, x2, 12), Rational(1, 2), 2, 12), 12))
D11inv = simpler(inverse(D11, 12))
print_var("D11", D11)
print_var("D11inv", D11inv)
# This right now is in fact D11 / (x1 - x0) = sqrt(1 - exp(-(x1 - x0)**2)) / (x1 - x0)
# Similarly, D11inv is in fact (x1 - x0) / sqrt(1 - exp(-(x1 - x0)**2))

# D20 = simpler(exponentiate(scalar_multiply(multiply(subtract(x0, x2, 8), subtract(x2, x2, 8), 8), Integer(-1)), 100000, 8))
# D20 = simpler(multiply(D20, D00inv, 8))
D20 = D00inv
print_var("D20", D20)

D21 = exp_series(11)
D21 = subtract(D21, one_series(11), 11)
D21.pop(0)
D21 = power_series(D21, scalar_multiply(multiply(subtract(x1, x0, 10), subtract(x2, x0, 10), 10), Integer(-1)), 2, 10)
D21 = multiply(D11inv, D21, 10)
D21 = multiply2(D21, 0, subtract(x2, x0, 11), 1, 11)
# D21 = multiply(D21, exponentiate(multiply(subtract(x1, x2, 7), subtract(x2, x2, 7), 7), 100000, 7), 7)
D21 = simpler(D21)
print_var("D21", D21)

D22squared = simpler(subtract(one_series(12), add(power(D20, 2, 0, 12), power(D21, 2, 1, 12), 12), 12))
# Simply: D22squared = 1 - D20**2 - D21**2
print_var("D22squared", D22squared)

dmean = add(scalar_multiply(D20, z0), scalar_multiply(D21, z1), 9)
dmean.insert(0, Integer(0))
dmean = simpler(subtract(x2, dmean, 10))
print_var("dmean", dmean)

D22squared.insert(0, Integer(0))
D22squared.insert(0, Integer(0))
error = add(power(dmean, 2, 4, 14), D22squared, 14)
error = simpler(error)
print_var("error", error) 

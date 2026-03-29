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

def substitute(terms, x, y):
	return [term.subs(x, y) for term in terms]

# ====================  SIMPLIFICATION  =========================

# We require MAX >= order_of_result - (pow - 1) * order_of_term
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

# We require MAX1 >= order_of_result // order_of_term
# In general, we require MAX2 >= order_of_result
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
# we require order_of_term >= 1
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
# We require MAX1 >= order_of_result + order2
# We require MAX2 >= order_of_result + 2*order2 - order1
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
# We require MAX >= order_of_result - order_of_term
def exponentiate_squared(terms, c, order_of_term, order_of_result):
	return power_series(exp2_series(c, order_of_result), terms, order_of_term, order_of_result)


# We require MAX1 >= order_of_result - order_of_term2
# We require MAX2 >= order_of_result - order_of_term1
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
# We require MAX >= order_of_terms//2 + order_of_result
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

ord = 6

a11, a12, a13, a14, a15, a16, a17, a18, a19, a110 = symbols('a11 a12 a13 a14 a15 a16 a17 a18 a19 a110')
z0 = symbols('z0')

a11 = z0
a12 = 0
a13 = -z0**3/2 - z0
a14 = 0


x0 = [Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)]
x1 = [Integer(0), a11, a12, a13, a14, a15, a16, a17, Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)]


# sum_k Cik * Cjk = exp(-|xi - xj|^2 / 2)
# Dik = Cik * exp(|xi - xN|^2 / 2)
# DNk = CNk
# sum_k Dik * Djk = exp(<xi - xN, xj - xN>)
# Aij = exp(<xi - xN, xj - xN>)


# D00 = exponentiate_squared(subtract(x0, x3, 9+ord), Rational(1, 2), 1, 10+ord)
D00inv = exponentiate_squared(subtract(x0, x1, 3+ord), Rational(-1, 2), 1, 4+ord)
D00inv = simpler(D00inv)
print_var("D00inv", D00inv)

D10 = D00inv
print_var("D10", D10)

D11squared = subtract(one_series(4+ord), multiply(D10, 0, D10, 0, 4+ord), 4+ord)
D11squared = simpler(D11squared)
print_var("D11squared", D11squared)



dmean = scalar_multiply(D10, z0)
dmean.insert(0, Integer(0))
dmean = simpler(subtract(x1, dmean, 2+ord))
print_var("dmean", dmean)

D11squared.insert(0, Integer(0))
D11squared.insert(0, Integer(0))
error = add(power(dmean, 2, 2, 4+ord), D11squared, 4+ord)
error = simpler(error)
print_var("error", error) 

print(str(error[4+ord]))


expr2 = diff(error[4+ord], a15)

sol = solve(expr2, a15)[0]
sol = simplify(sol)
pprint(sol)
pprint(str(sol))

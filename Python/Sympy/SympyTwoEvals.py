from sympy import *
init_printing()

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

def one_plus_series(term, order):
	result = [Integer(1)]
	result.append(term)
	for k in range(2, order + 1):
		result.append(Integer(0))
	return result

def zero_plus_series(term, order):
	result = [Integer(0)]
	result.append(term)
	for k in range(2, order + 1):
		result.append(Integer(0))
	return result

def scalar_series(term, order):
	result = [term]
	for k in range(1, order + 1):
		result.append(Integer(0))
	return result

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


# we assume that terms[0] = 1
def inverse(terms, order_of_result):
	result = [Integer(0) for k in range(order_of_result + 1)]
	result[0] = Integer(1)
	for s in range(order_of_result + 1):
		for i in range(s):
			result[s] -= terms[s - i] * result[i]
	return result


# we assume that terms[0] = 1
def sqrt(terms, order_of_result):
	result = [Integer(0) for k in range(order_of_result + 1)]
	result[0] = Integer(1)
	if order_of_result > 0:
		result[1] = Rational(1, 2) * terms[1]
	for s in range(2, order_of_result + 1):
		result[s] = terms[s]
		for i in range(1, s):
			result[s] -= result[s - i] * result[i]
		result[s] *= Rational(1, 2)
	return result

# =====================================================================================================

# We expand stuff in sqrt(p).
order = 4
mu = symbols('mu')
# K11, K12, K13, K22, K23, K33 = symbols('K11 K12 K13 K22 K23 K33')
# E12, E13, E23 = symbols('E12 E13 E23')
# X211, X212, X213, X214 = symbols('X211 X212 X213 X214')
# X311, X312, X313, X314 = symbols('X311 X312 X313 X314')
K11, K12, K13, K22, K33 = symbols('K11 K12 K13 K22 K33')
E23 = symbols('E23')
X212, X213, X214, X215, X216 = symbols('X212 X213 X214 X215 X216')
X313, X314, X315, X316 = symbols('X313 X314 X315 X316')

# Results:
X311 = K13 / mu
X312 = Integer(0)

# Possibly flip the sign in BOTH of the below!
K23 = -K13
X211 = -E23 / mu

X21 = [Integer(0), X211, X212, X213, X214, X215, X216]
X31 = [Integer(0), X311, X312, X313, X314, X315, X316]

exp_X21 = exponentiate_squared(X21, Rational(-1, 2), 1, order)
exp_X31 = exponentiate_squared(X31, Rational(-1, 2), 1, order)
exp_diff = exponentiate_squared(subtract(X21, X31, order), Rational(-1, 2), 1, order)

# We calculate the optimum with assumption that both evals occur at 1.
sqrt_p_C11 = sqrt(add(one_series(order), [Integer(0), Integer(0)] + scalar_series(K11, order - 2), order), order)
sqrt_inv_p_C21 = divide(scalar_multiply(exp_X21, K12), 0, sqrt_p_C11, 0, order)
sqrt_inv_p_C31_without_X32 = divide(scalar_multiply(exp_X31, K13), 0, sqrt_p_C11, 0, order)

p_K22_prim = add(one_series(order), [Integer(0), Integer(0)] + subtract(one_series(order - 2), scalar_multiply(multiply(exp_X21, 0, exp_X21, 0, order - 2), Integer(1) - K22), order - 2), order)
sqrt_p_C22 = sqrt(subtract(p_K22_prim, [Integer(0), Integer(0), Integer(0), Integer(0)] + multiply(sqrt_inv_p_C21, 0, sqrt_inv_p_C21, 0, order - 4), order), order)
K23prim_without_X32 = subtract(scalar_multiply(exp_diff, E23), scalar_multiply(multiply(exp_X21, 0, exp_X31, 0, order), E23 - K23), order)
sqrt_inv_p_C32_without_X32 = divide(subtract(K23prim_without_X32, [Integer(0), Integer(0)] + multiply(sqrt_inv_p_C21, 0, sqrt_inv_p_C31_without_X32, 0, order - 2), order), 0, sqrt_p_C22, 0, order)

sqrt_inv_p_q = divide(scalar_multiply(sqrt_inv_p_C32_without_X32, Integer(2)), 0, add(scalar_multiply(exp_X31, mu), [Integer(0), Integer(0)] + multiply(X31[1:], 0, sqrt_inv_p_C31_without_X32, 0, order - 2), order), 0, order)
sqrt_inv_p_t = divide(subtract(sqrt(add(one_series(order + 2), [Integer(0), Integer(0)] + multiply(sqrt_inv_p_q, 0, sqrt_inv_p_q, 0, order), order + 2), order + 2), one_series(order + 2), order + 2)[2:], 0, sqrt_inv_p_q, 0, order)

exp_t = exponentiate_squared([Integer(0)] + sqrt_inv_p_t, Rational(-1, 2), 1, order)
result = multiply(exp_t, 0, add(scalar_multiply(exp_X31, mu), [Integer(0), Integer(0)] + add(multiply(X31[1:], 0, sqrt_inv_p_C31_without_X32, 0, order - 2), multiply(sqrt_inv_p_C32_without_X32, 0, sqrt_inv_p_t, 0, order - 2), order - 2), order), 0, order)

pprint(simpler(result))
# pprint(sqrt_p_C11)
# pprint(simplify(result[order]))

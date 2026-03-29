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


# =============================================

order = 5

a2, a3, a4, a5, a6, a7 = symbols('a2 a3 a4 a5 a6 a7')
b2, b3, b4, b5, b6, b7 = symbols('b2 b3 b4 b5 b6 b7')
c = symbols('c')

f = [Integer(0), Integer(1), a2, a3, a4, a5, a6, a7]
# hinverse = [Integer(0), Integer(1), b2, b3, b4, b5, b6, b7]

# f_of_hinverse = simpler(power_series(f, hinverse, 1, order + 1))

# solution2 = simplify(b2 - f_of_hinverse[2])
# print_var("b2:", solution2)
# solution2 /= (1 - c)
# f_of_hinverse = substitute(f_of_hinverse, b2, solution2)

# solution3 = simplify(b3 - f_of_hinverse[3])
# print_var("b3:", solution3)
# solution3 /= (1 - c**2)
# f_of_hinverse = substitute(f_of_hinverse, b3, solution3)

# solution4 = simplify(b4 - f_of_hinverse[4])
# print_var("b4:", solution4)
# solution4 /= (1 - c**3)
# f_of_hinverse = substitute(f_of_hinverse, b4, solution4)

# print_var("b4 c:", simplify(solution4.subs(c, Integer(0))))

# solution5 = simplify(b5 - f_of_hinverse[5])
# print_var("b5:", solution5)
# solution5 /= (1 - c**4)
# f_of_hinverse = substitute(f_of_hinverse, b5, solution5)

# print_var("b5 c:", simplify(solution5.subs(c, Integer(0))))


# solution6 = simplify(b6 - f_of_hinverse[6])
# print_var("b6:", solution6)
# solution6 /= (1 - c**5)
# f_of_hinverse = substitute(f_of_hinverse, b6, solution6)

# print_var("b6 c:", simplify(solution6.subs(c, Integer(0))))




f2 = power_series(f, f, 1, order)
f3 = power_series(f2, f, 1, order)
f4 = power_series(f3, f, 1, order)
ones_series = [Integer(0) for i in range(order + 1)]
m_ones_series = [Integer(0) for i in range(order + 1)]
ones_series[1] = Integer(1)
m_ones_series[1] = Integer(-1)
# combination = add(f3, add(scalar_multiply(f2, Integer(-3)), add(scalar_multiply(f, Integer(3)), m_ones_series, order), order), order)
combination = add(f4, add(scalar_multiply(f3, Integer(-4)), add(scalar_multiply(f2, Integer(6)), add(scalar_multiply(f, Integer(-4)), ones_series, order), order), order), order)
combination = simpler(combination)

print_var("0", combination[0])
print_var("1", combination[1])
print_var("2", combination[2])
print_var("3", combination[3])
print_var("4", combination[4])
print_var("5", combination[5])

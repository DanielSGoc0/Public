from sympy import *
from sympy.core.sympify import sympify
import random


def to_rational(x):
	# sympify handles ints, Fractions, "a/b", etc.
	x = sympify(x)
	# Force exact rational when possible (e.g. 0.5 -> 1/2).
	# If you pass floats like 0.1, this will become 1/10 only if written as Rational("0.1")
	return Rational(x)

def in_span_rational(vectors, w):
    """
    Exact check whether w ∈ span_Q{v1,...,vn} (i.e., rational coefficients).

    Inputs may be ints, Fractions, strings like "3/7", or sympy Rationals.
    Returns (ok, coeffs) where coeffs is one exact solution in Q^n if it exists.
    """

    w = Matrix([to_rational(x) for x in w])
    cols = [Matrix([to_rational(x) for x in v]) for v in vectors]
    V = Matrix.hstack(*cols)  # d x n

    # Solve V*c = w exactly. linsolve returns a set of solutions (possibly empty).
    V.gauss_jordan_solve(w)  # returns (particular_solution, param_matrix)

    return True


def in_span_rational_safe(vectors, w):
    """Same as above, but returns (False, None) instead of throwing if inconsistent."""
    try:
        ok = in_span_rational(vectors, w)
        return ok
    except Exception:
        return False


# --- Example ---
if __name__ == "__main__":
	for attempt in range(10000):

		v = [sympify(random.randint(-2, 2)) for i in range(5)]
		v[4] = sympify("1") - (v[0] + v[1] + v[2] + v[3])

		a = sympify(random.randint(-2, 2))
		b = sympify(random.randint(-2, 2))
		c = None
		d = None
		e = None

		M = Matrix([[v[2], v[3], v[4]], [v[4], v[0], v[1]], [v[0], v[1], v[2]]])
		u = Matrix([-v[0]*b - v[1]*a,  -v[2]*b - v[3]*a, -v[3]*b - v[4]*a])

		sol_set = linsolve((M, u))
		if not sol_set:
			# print("No solution")
			continue
		else:
			sol_tuple = next(iter(sol_set))

			params = set().union(*(expr.free_symbols for expr in sol_tuple))
			if not params:
				# unique solution
				(c, d, e) = sol_tuple
			else:
				subs_dict = {}
				for p in params:
					num = random.randint(-3, 3)
					subs_dict[p] = sympify(num)
				
				(c, d, e) = [expr.subs(subs_dict) for expr in sol_tuple]

				print(a, b, c, d, e)
		

		# print(v)
		# print(a, b, c, d, e)


		# M1 = [[a[2], a[3], a[4], "1"], [a[3], a[4], a[0], "1"], [a[4], a[0], a[1], "1"], [a[0], a[1], a[2], "1"], [a[1], a[2], a[3], "1"]]
		# M2 = [[a[1], a[3], a[4], "1"], [a[2], a[4], a[0], "1"], [a[3], a[0], a[1], "1"], [a[4], a[1], a[2], "1"], [a[0], a[2], a[3], "1"]]
            
		M1 = [[b, d, e, "1"], [a, c, d, "1"], [e, b, c, "1"], [d, a, b, "1"], [c, e, a, "1"]]
		M2 = [[c, d, e, "1"], [b, c, d, "1"], [a, b, c, "1"], [e, a, b, "1"], [d, e, a, "1"]]

		w = ["0", "0", "0", "1"]

		ok1 = in_span_rational_safe(M1, w)
		ok2 = in_span_rational_safe(M2, w)
		if ok1 != ok2:
			print("DONE")
			print(a)
			break
		if not ok1:
			print(v)
			print(a, b, c, d, e)
          
		print(attempt, ok1, ok2)

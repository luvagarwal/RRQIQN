import random
from sympy import symbols, conjugate, Abs, combsimp, latex


def diagonalElementsInLatex(level):
    def sharelatex(y):
        y = y.replace('\\left(', '(')
        y = y.replace('\\right)', ')')
        return y

    m = print_density_matrix(level)
    for i in xrange(4):
        print sharelatex(latex(m[i, i])), "\n\n"

def get_resource_parameters(level=None):
    if level is not None:
        return symbols('p%s l%s L%s P%s' % tuple([level + 1] * 4), real=True, positive=True)
    return symbols('p l L P', real=True, positive=True)

def get_basis_parameters(level=None):
    if level is not None:
        return symbols('q%s Q%s' % tuple([level + 1] * 2), real=True, positive=True)
    return symbols('q Q', real=True, positive=True)

def get_initial_secret_parameters():
    return symbols('a b', real=True, positive=True)

def subsitute(expr):
    p, l, L, P = get_resource_parameters()
    q, Q = get_basis_parameters()
    a, b = get_initial_secret_parameters()
    change = [l, a, p, q, b]
    for sym in change:
        expr = expr.subs(sym*conjugate(sym), Abs(sym)**2)
    return expr

def simplify_density_matrix(density_matrix):
    for row in xrange(4):
        for col in xrange(4):
            density_matrix[row, col] = subsitute(combsimp(density_matrix[row, col]))
    return density_matrix

def debug(**kwargs):
    def upto_n_decimal_points(num, upto=2):
        num = str(num).split('.')
        return "%s.%s" % (num[0], num[1][:2])

    for key, val in kwargs.iteritems():
        print "%s = %s" % (key, upto_n_decimal_points(val))

def monte_carlo(fn, var, domain=(0, 1), rounds=500):
    total = 0
    start, end = domain
    for _ in xrange(rounds):
        random_val = start + random.random() * (end - start)
        total += fn.subs(var, random_val)
    return total / rounds

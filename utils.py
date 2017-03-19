from sympy import *
from sympy.matrices import *


def diagonalElementsInLatex(level):
    def sharelatex(y):
        y = y.replace('\\left(', '(')
        y = y.replace('\\right)', ')')
        return y

    m = print_density_matrix(level)
    for i in xrange(4):
        print sharelatex(latex(m[i, i]))
        print
        print

# def compact(a):
#     return collect(a, a.atoms(Symbol), func=factor)

def subsitute2(expr):
    p, l, q, Q, L, P, b, a = symbols('p l q Q L P, b, a')
    change = [l, a, p, q, b]
    for sym in change:
        expr = expr.subs(sym*conjugate(sym), Abs(sym)**2)
    return expr

def subsitute(coeff_comm):
    p, l, q, Q, L, P = symbols('p l q Q L P')
    coeff_comm = coeff_comm.subs(L, 1/(1+Abs(l)**2)**0.5).subs(P, 1/(1+Abs(p)**2)**0.5).subs(Q, 1/(1+Abs(q)**2)**0.5)
    return coeff_comm

def subsitute3(expr):
    p, l, q, Q, L, P, b, a = symbols('p l q Q L P, b, a')
    change = [l, a, p, q, b]
    for sym in change:
        expr = expr.subs()

def simplify_density_matrix(density_matrix):
    for row in xrange(4):
        for col in xrange(4):
            density_matrix[row, col] = subsitute2(combsimp(density_matrix[row, col]))
    return density_matrix

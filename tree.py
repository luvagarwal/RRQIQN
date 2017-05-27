from __future__ import division
import re
from argparse import Namespace
from sympy import pprint
from sympy.matrices import Matrix
from sympy.parsing.sympy_parser import parse_expr
from sympy import sqrt

from utils import *


TYPES = {
    'PHI': 'PHI',
    'PSI': 'PSI',
}
TYPES = Namespace(**TYPES)


def state_factory(type, level, coeffs, coeff_comm=1):
    if type == TYPES.PHI:
        return StatePhi(level, coeffs, coeff_comm)
    elif type == TYPES.PSI:
        return StatePsi(level, coeffs, coeff_comm)


class State(object):
    def __init__(self, level, coeffs, coeff_comm=1):
        self.coeffs = coeffs
        self.coeff_comm = coeff_comm
        self.resource_parameters = get_resource_parameters()
        self.basis_parameters = get_basis_parameters()
        self.initial_secret_params = get_initial_secret_parameters()
        self.level = level

    def get_children(self, fall_type=None):
        try:
            self.children
        except AttributeError:
            self.children = self.lazy_load_children()

        if fall_type is None:
            return self.children

        mapping = {
            'phi_plus': 0,
            'phi_minus': 1,
            'psi_plus': 2,
            'psi_minus': 3,
        }
        return self.children[mapping[fall_type]]

    def get_descendants(self, depth):
        if depth < 0:
            raise Exception("Depth can't be negative")
        if depth == 0:
            return [self]

        descendants = []
        for child in self.get_children():
            descendants += child.get_descendants(depth - 1)
        return descendants

    def print_descendants(self, depth):
        descendants = self.get_descendants(depth)
        for state in descendants:
            state.pprint()

    def as_density_matrix(self):
        " Represent this pure state as density matrix "
        coeffs = self.get_coefficients_as_general_state()
        matrix = Matrix(coeffs)
        matrix = matrix * matrix.conjugate().transpose()
        coeff_comm = self.coeff_comm
        matrix = matrix * coeff_comm * coeff_comm
        return matrix

    def get_density_matrix(self, depth):
        nodes = self.get_descendants(depth)
        density_matrix = Matrix([[0]*4 for i in xrange(4)])

        for node in nodes:
            matrix = node.as_density_matrix()
            density_matrix += matrix

        return simplify_density_matrix(density_matrix)

    def print_density_matrix(self, depth, only_diagonal=True):
        density_matrix = self.get_density_matrix(depth)
        for row in xrange(4):
            for col in xrange(4):
                if only_diagonal and row != col:
                    continue
                pprint(density_matrix[row, col])

    def as_single_qubit(self):
        a_dash, b_dash = self.coeffs
        norm = sqrt(a_dash**2 + b_dash**2)
        a_dash /= norm
        b_dash /= norm
        p, l, _, _ = self.resource_parameters
        q, _ = self.basis_parameters
        return SingleQubit([a_dash, b_dash], self.initial_secret_params, [p, l, q])

    def get_overlap(self):
        pass

    def pprint(self):
        pprint(self.as_sympy_expr())

    def __repr__(self):
        return repr(self.as_sympy_expr())


class StatePhi(State):
    def lazy_load_children(self):
        a, b = self.coeffs
        c = self.coeff_comm
        level = self.level
        p, l, L, P = self.resource_parameters
        q, Q = self.basis_parameters
        return  [
            StatePhi(level + 1, [a, b*q*l.conjugate()], c*Q*L),
            StatePhi(level + 1, [a*l, -b*q], c*Q*L),
            StatePsi(level + 1, [a*q, b*p.conjugate()], c*Q*P),
            StatePsi(level + 1, [a*q*p, -b], c*Q*P)
        ]

    def as_sympy_expr(self):
        expr = "%s*(%s*Symbol('|00|') + %s*Symbol('|11|'))" % (self.coeff_comm, self.coeffs[0], self.coeffs[1])
        expr = parse_expr(expr)
        return expr

    def get_coefficients_as_general_state(self):
        return [self.coeffs[0], 0, 0, self.coeffs[1]]


class StatePsi(State):
    def lazy_load_children(self):
        a, b = self.coeffs
        c = self.coeff_comm
        level = self.level
        p, l, L, P = self.resource_parameters
        q, Q = self.basis_parameters
        return  [
            StatePsi(level + 1, [b*q*l.conjugate(), a], c*Q*L),
            StatePsi(level + 1, [-b*q, a*l], c*Q*L),
            StatePhi(level + 1, [b*p.conjugate(), a*q], c*Q*P),
            StatePhi(level + 1, [-b, a*q*p], c*Q*P)
        ]

    def as_sympy_expr(self):
        expr = "%s*(%s*Symbol('|01|') + %s*Symbol('|10|'))"%(self.coeff_comm, self.coeffs[0], self.coeffs[1])
        expr = parse_expr(expr)
        return expr

    def get_coefficients_as_general_state(self):
        return [0, self.coeffs[0], self.coeffs[1], 0]


class SingleQubit(object):
    def __init__(self, coeffs, initial_secret_params, resouce_and_basis_params):
        self.coeffs = coeffs
        self.initial_secret_params = initial_secret_params
        self.resouce_and_basis_params = resouce_and_basis_params

    def subs_param(self, param, val):
        self.coeffs = [coeff.subs(param, val) for coeff in self.coeffs]

    def subs_params(self, params, with_val=False):
        if not with_val:
            random_params = {}
            for param in params:
                random_params[param] = random.random()
        else:
            random_params = params

        for param, val in random_params.items():
            self.subs_param(param, val)
        return random_params

    def subs_all(self):
        a, b = self.initial_secret_params
        self.subs_b_by_a()
        return self.subs_params(self.resouce_and_basis_params + [a])

    def subs_b_by_a(self):
        a, b = self.initial_secret_params
        self.subs_param(b, sqrt(1 - a**2))

    def subs_secret_params(self):
        a, b = self.initial_secret_params
        self.subs_b_by_a()
        return self.subs_params([a])

    def subs_resource_and_basis_params(self):
        return self.subs_params(self.resouce_and_basis_params)

    def get_abs_coeffs(self):
        abs_coeff = lambda coeff: -1*coeff if -1 in coeff.args else coeff
        return map(abs_coeff, self.coeffs)

    def apply_basic_unitary(self):
        a, b = self.initial_secret_params
        coeffs = self.get_abs_coeffs()
        coeffs = coeffs if a in coeffs[0].args or a is coeffs[0] else coeffs[::-1]
        return SingleQubit(coeffs, self.initial_secret_params, self.resouce_and_basis_params)

    def apply_advanced_unitary(self):
        a, b = self.initial_secret_params
        self.subs_b_by_a()
        b = (1 - a**2) ** 0.5
        a_dash, b_dash = self.coeffs

        alpha = a*a_dash + b*b_dash
        avg_alpha = monte_carlo(alpha, a)
        debug(avg_alpha=avg_alpha)

        beta = a_dash*b - b_dash*a
        avg_beta = monte_carlo(beta, a)
        debug(avg_beta=avg_beta)

        x = avg_alpha / (avg_alpha**2 + avg_beta**2)**0.5
        debug(x=x)
        y = (1 - x**2) ** 0.5
        optimized_coeffs = [x * a_dash - y * b_dash, y * a_dash + x * b_dash]
        return SingleQubit(optimized_coeffs, self.initial_secret_params, self.resouce_and_basis_params)

    def norm(self):
        a, b = self.coeffs
        return (a**2 + b**2) ** 0.5

    def dot_product(self, state):
        a1, b1 = self.coeffs
        a2, b2 = state.coeffs
        return (a1*a2 + b1*b2) / (self.norm() * state.norm())

    def as_sympy_expr(self):
        expr = "%s*Symbol('|0|') + %s*Symbol('|1|')"%(self.coeffs[0], self.coeffs[1])
        expr = parse_expr(expr)
        return expr

    def pprint(self):
        pprint(self.as_sympy_expr())



def cache(fn):
    pass

def create_root():
    " get root node "
    a, b = get_initial_secret_parameters()
    return state_factory(TYPES.PHI, 0, [a, b])


def main():
    pass

if __name__ == "__main__":
    root = create_root()
    root.print_density_matrix(2)

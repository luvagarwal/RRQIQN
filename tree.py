from __future__ import division
import re
from argparse import Namespace
from sympy import *
from sympy.matrices import *
from sympy.parsing.sympy_parser import parse_expr

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
        self.level = level
        self.resource_parameters = get_resource_parameters()
        self.basis_parameters = get_basis_parameters()

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
        density_matrix = Matrix([[0, 0, 0, 0] for i in xrange(4)])

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

    def get_abs_coeffs(self):
        abs_coeff = lambda coeff: -1*coeff if -1 in coeff.args else coeff
        return map(abs_coeff, self.coeffs)

    def as_single_qubit(self):
        """
        1. Convert 2 qubit to 1 qubit
        2. Apply basic unitary operator
        """
        a, b = get_initial_secret_parameters()
        coeffs = self.get_abs_coeffs()
        return coeffs if a in coeffs[0].args or a is coeffs[0] else coeffs[::-1]

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

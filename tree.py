#!/usr/bin/env python
from __future__ import division
import re
from argparse import Namespace
from sympy import *
from sympy.matrices import *
from sympy.physics.quantum import TensorProduct
from sympy.parsing.sympy_parser import parse_expr

from utils import *


TYPES = {
    'PHI': 'PHI',
    'PSI': 'PSI',
}
TYPES = Namespace(**TYPES)


def state_factory(type, level, coeff, coeff_comm=1):
    if type == TYPES.PHI:
        return StatePhi(level, coeff, coeff_comm)
    elif type == TYPES.PSI:
        return StatePsi(level, coeff, coeff_comm)


class State(object):
    def __init__(self, level, coeff, coeff_comm=1):
        self.coeff = coeff
        self.coeff_comm = coeff_comm
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
            raise Exception('Depth can\'t be negative')
        if depth == 0:
            return [self]

        descendants = []
        for child in self.get_children():
            descendants += child.get_descendants(depth - 1)
        return descendants

    def as_density_matrix(self):
        " Represent this pure state as density matrix "
        coeffs = self.get_coefficients_as_general_state()
        matrix = Matrix(coeffs)
        matrix = matrix * matrix.conjugate().transpose()
        coeff_comm = self.coeff_comm
        matrix = matrix * coeff_comm * coeff_comm
        return matrix

    def __repr__(self):
        return repr(self.as_sympy_expr())


class StatePhi(State):
    def lazy_load_children(self):
        a, b = self.coeff
        c = self.coeff_comm
        level = self.level
        q, Q = symbols('q%s Q%s' % (level+1, level+1))
        p, l, L, P = symbols('p l L P')
        return  [
            StatePhi(level + 1, [a, b*q*l.conjugate()], c*Q*L),
            StatePhi(level + 1, [a*l, -b*q], c*Q*L),
            StatePsi(level + 1, [a*q, b*p.conjugate()], c*Q*P),
            StatePsi(level + 1, [a*q*p, -b], c*Q*P)
        ]

    def as_sympy_expr(self):
        expr = "%s*(%s*Symbol('|00|') + %s*Symbol('|11|'))" % (self.coeff_comm, self.coeff[0], self.coeff[1])
        expr = parse_expr(expr)
        return expr

    def get_coefficients_as_general_state(self):
        return [self.coeff[0], 0, 0, self.coeff[1]]


class StatePsi(State):
    def lazy_load_children(self):
        a, b = self.coeff
        c = self.coeff_comm
        level = self.level
        q, Q = symbols('q%s Q%s' % (level+1, level+1))
        p, l, L, P = symbols('p l L P')
        return  [
            StatePsi(level + 1, [b*q*l.conjugate(), a], c*Q*L),
            StatePsi(level + 1, [-b*q, a*l], c*Q*L),
            StatePhi(level + 1, [b*p.conjugate(), a*q], c*Q*P),
            StatePhi(level + 1, [-b, a*q*p], c*Q*P)
        ]

    def as_sympy_expr(self):
        expr = "%s*(%s*Symbol('|01|') + %s*Symbol('|10|'))"%(self.coeff_comm, self.coeff[0], self.coeff[1])
        expr = parse_expr(expr)
        return expr

    def get_coefficients_as_general_state(self):
        return [0, self.coeff[0], self.coeff[1], 0]


def create_root():
    " get root node "
    a, b = symbols('a b')
    return state_factory(TYPES.PHI, 0, [a, b])

root = create_root()


def get_density_matrix(level):
    nodes = root.get_descendants(level)
    density_matrix = Matrix([[0, 0, 0, 0] for i in xrange(4)])

    for node in nodes:
        matrix = node.as_density_matrix()
        density_matrix += matrix

    return simplify_density_matrix(density_matrix)

def print_density_matrix(level):
    density_matrix = get_density_matrix(level)
    for row in xrange(4):
        for col in xrange(4):
            pprint(density_matrix[row, col])

def get_states(level):
    states = root.get_descendants(level)
    states = [state.as_sympy_expr() for state in states]
    return states

def print_states(level):
    states = get_states(level)
    for state in states:
        pprint(state)

def main():
    pass

if __name__ == "__main__":
    print get_states(5)
    print print_density_matrix(5)

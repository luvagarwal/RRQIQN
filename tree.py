from __future__ import division
import re
from argparse import Namespace
from sympy import pprint
from sympy.matrices import Matrix
from sympy.parsing.sympy_parser import parse_expr
from sympy import sqrt

from utils import *
from paramvalues import ParamsSingleton


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
    """
    Representation of a general 2 qubit state. This class is
    inherited by classes corresponding to specific 2 qubit states
    like PHI and PSI state.
    """

    def __init__(self, level, coeffs, coeff_comm=1):
        # It's an array of size 2 with a_dash and b_dash as values
        self.coeffs = coeffs
        # It is the square root of probablity of this state's occurence.
        self.coeff_comm = coeff_comm
        # p, l
        self.resource_parameters = get_resource_parameters()
        # q
        self.basis_parameters = get_basis_parameters()
        # a, b
        self.initial_secret_params = get_initial_secret_parameters()
        # number of swappings
        self.level = level

    def get_children(self, fall_type=None):
        """
        Return it's 4 children states after a swapping, arranged
        according to the measurement fall as:
        [phi_plus, phi_minus, psi_plus, psi_minus]

        Finding children nodes is an expensive process specially if
        we are dealing with too many swappings. So it will only
        calculate children on demand and cache them. It doesn't
        precalculate them.
        """
        try:
            self.children
        except AttributeError:
            self.children = self.lazy_load_children()

        return self.children

    def get_children_with_fall_type(fall_type):
        """
        Return children which is generated when measurement falls
        to `fall_type`.
        """
        mapping = {
            'phi_plus': 0,
            'phi_minus': 1,
            'psi_plus': 2,
            'psi_minus': 3,
        }

        return self.get_children()[mapping[fall_type]]

    def get_descendants(self, depth):
        """
        Return states after swappings = `depth`. So it will
        return array of size 4^(depth)
        """
        if depth < 0:
            raise Exception("Depth can't be negative")
        if depth == 0:
            return [self]

        descendants = []
        for child in self.get_children():
            descendants += child.get_descendants(depth - 1)
        return descendants

    def print_descendants(self, depth):
        " Pretty states after no of swappings = `depth` "
        descendants = self.get_descendants(depth)
        for state in descendants:
            state.pprint()

    def as_density_matrix(self):
        """
        NOT USED ANYMORE (but works fine)
        Represent this pure state as density matrix.
        """
        coeffs = self.get_coefficients_as_general_state()
        matrix = Matrix(coeffs)
        matrix = matrix * matrix.conjugate().transpose()
        coeff_comm = self.coeff_comm
        matrix = matrix * coeff_comm * coeff_comm
        return matrix

    def get_density_matrix(self, depth):
        """
        Represent all the states after swappings = `depth`
        as density matrix.
        """
        nodes = self.get_descendants(depth)
        density_matrix = Matrix([[0]*4 for i in xrange(4)])

        for node in nodes:
            matrix = node.as_density_matrix()
            density_matrix += matrix

        return simplify_density_matrix(density_matrix)

    def print_density_matrix(self, depth, only_diagonal=True):
        """
        Pretty print density matrix after swappings = `depth`.
        Since non diagonal elements are zero, it prints only diagonal
        elements by default.
        """
        density_matrix = self.get_density_matrix(depth)
        for row in xrange(4):
            for col in xrange(4):
                if only_diagonal and row != col:
                    continue
                pprint(density_matrix[row, col])

    def as_single_qubit(self):
        """
        Reconstruct this state to single qubit and return that.
        """
        a_dash, b_dash = self.coeffs
        norm = sqrt(a_dash**2 + b_dash**2)
        a_dash /= norm
        b_dash /= norm
        p, l, _, _ = self.resource_parameters
        q, _ = self.basis_parameters
        return SingleQubit([a_dash, b_dash], self.initial_secret_params, [p, l, q])

    def pprint(self):
        pprint(self.as_sympy_expr())

    def __repr__(self):
        return repr(self.as_sympy_expr())


class StatePhi(State):
    """
    It includes specific functions corresponding to nodes in PHI state i.e. (|00\, |11\)
    representation. In general, a `State` can either be in  (|00\, |11\) or (|01\, |10\)
    representation. This class contains functions corresponding to (|00\, |11\)
    representation.
    """
    def lazy_load_children(self):
        " Calculate it's children after a swapping "
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
        """
        Used to print this node in a pretty format. It is responsible for those beautiful
        states that you see on the console.
        """
        expr = "%s*(%s*Symbol('|00|') + %s*Symbol('|11|'))" % (self.coeff_comm, self.coeffs[0], self.coeffs[1])
        expr = parse_expr(expr)
        return expr

    def get_coefficients_as_general_state(self):
        """
        Return coefficients assuming this state in
        x|00\ + y|01\ + z|10\ + w|11\ . Since this is a PHI state,
        it will have y and z as zero.
        """
        return [self.coeffs[0], 0, 0, self.coeffs[1]]


class StatePsi(State):
    """
    It includes specific functions corresponding to nodes in PSI state i.e. (|01\, |10\)
    representation. In general, a `State` can either be in  (|00\, |11\) or (|01\, |10\)
    representation. This class contains functions corresponding to (|01\, |10\)
    representation.
    """
    def lazy_load_children(self):
        " Calculate it's children after a swapping "
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
        """
        Used to print this node in a pretty format. It is responsible for those beautiful
        states that you see on the console.
        """
        expr = "%s*(%s*Symbol('|01|') + %s*Symbol('|10|'))"%(self.coeff_comm, self.coeffs[0], self.coeffs[1])
        expr = parse_expr(expr)
        return expr

    def get_coefficients_as_general_state(self):
        """
        Return coefficients assuming this state in
        x|00\ + y|01\ + z|10\ + w|11\ . Since this is a PSI state,
        it will have x and w as zero.
        """
        return [0, self.coeffs[0], self.coeffs[1], 0]




class SingleQubit(object):
    """
    This class represents a single qubit state.
    """
    def __init__(self, coeffs, initial_secret_params, resource_and_basis_parameters):
        # [a_dash, b_dash] for (a_dash|0\ + b_dash|1\). Note that for initial
        # secret a_dash = a and b_dash = b
        self.coeffs = coeffs
        # [a, b] (initial secret params)
        self.initial_secret_params = initial_secret_params
        # [p, l, q]
        self.resource_and_basis_parameters = resource_and_basis_parameters

    def subs_param(self, param, val):
        " Substitute `param` = `val` in it's coefficients "
        self.coeffs = [coeff.subs(param, val) for coeff in self.coeffs]

    def subs_params(self, params, with_val=False):
        """
        If with_val = False (by default) then:
        - `params` is a list of parameters
        - Randomly substitute these params. For eg. if we want to substitute p and q randomly,
         `params` will be [p, q].
        If with_val = True:
        - `params` is a dictionary with params as keys and values as their substitution values
        """
        if not with_val:
            random_params = ParamsSingleton().get_params_values(params)
        else:
            random_params = params

        for param, val in random_params.items():
            self.subs_param(param, val)
        return random_params

    def subs_all(self):
        """
        Substitute all parameters i.e. [a, b, p, q, l]. Note that b is substituted with
        (1-a**2)**0.5
        """
        a, b = self.initial_secret_params
        self.subs_b_by_a()
        return self.subs_params(self.resource_and_basis_parameters + [a])

    def subs_b_by_a(self):
        " Substitute b by (1-a**2)**0.5 "
        a, b = self.initial_secret_params
        self.subs_param(b, sqrt(1 - a**2))

    def subs_secret_params(self):
        " Substitute a randomly and b accordingly "
        a, b = self.initial_secret_params
        self.subs_b_by_a()
        return self.subs_params([a])

    def subs_resource_and_basis_params(self):
        " Substitute p, q, l randomly "
        return self.subs_params(self.resource_and_basis_parameters)

    def get_abs_coeffs(self):
        abs_coeff = lambda coeff: -1*coeff if -1 in coeff.args else coeff
        return map(abs_coeff, self.coeffs)

    def apply_basic_unitary(self):
        """
        Convert this state to standard a|0\ + b|1\ state like the way it is done
        for maximally entangled state after reconstruction.
        """
        a, b = self.initial_secret_params
        coeffs = self.get_abs_coeffs()
        coeffs = coeffs if a in coeffs[0].args or a is coeffs[0] else coeffs[::-1]
        return SingleQubit(coeffs, self.initial_secret_params, self.resource_and_basis_parameters)

    def calculate_advanced_unitary(self):
        a, b = self.initial_secret_params
        self.subs_b_by_a()
        b = (1 - a**2) ** 0.5
        a_dash, b_dash = self.coeffs

        alpha = a*a_dash + b*b_dash
        # avg_alpha = monte_carlo(alpha, a)
        # debug(avg_alpha=avg_alpha)

        beta = a_dash*b - b_dash*a
        # avg_beta = monte_carlo(beta, a)
        # debug(avg_beta=avg_beta)

        x = monte_carlo(alpha / (alpha**2 + beta**2) ** 0.5, a)
        # x = avg_alpha / (avg_alpha**2 + avg_beta**2)**0.5
        # debug(x=x)
        y = (1 - x**2) ** 0.5
        return x, y

    def apply_advanced_unitary(self, return_unitary=False):
        " Apply advanced unitary "
        x, y = self.calculate_advanced_unitary()
        a_dash, b_dash = self.coeffs
        optimized_coeffs = [x * a_dash - y * b_dash, y * a_dash + x * b_dash]
        new_state = SingleQubit(optimized_coeffs, self.initial_secret_params, self.resource_and_basis_parameters)
        if return_unitary:
            return x, new_state
        return new_state

    def norm(self):
        " Return it's norm which is square root of sum of it's squared coefficients "
        a, b = self.coeffs
        return (a**2 + b**2) ** 0.5

    def dot_product(self, state):
        " Return dot product of itself with the given `state` "
        a1, b1 = self.coeffs
        a2, b2 = state.coeffs
        return (a1*a2 + b1*b2) / (self.norm() * state.norm())

    def as_sympy_expr(self):
        """
        Used to print this node in a pretty format. It is responsible for those beautiful
        states that you see on the console.
        """
        expr = "%s*Symbol('|0|') + %s*Symbol('|1|')"%(self.coeffs[0], self.coeffs[1])
        return parse_expr(expr)

    def __repr__(self):
        return repr(self.as_sympy_expr())

    def pprint(self):
        pprint(self.as_sympy_expr())


def create_root():
    """
    Create root node
    Assumption: we start with a|00\ + b|11\ .
    """
    a, b = get_initial_secret_parameters()
    return state_factory(TYPES.PHI, 0, [a, b])


def main():
    pass

if __name__ == "__main__":
    root = create_root()
    root.print_descendants(2)

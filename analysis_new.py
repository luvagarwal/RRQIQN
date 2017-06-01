from __future__ import division
import random
from sympy import sqrt

from utils import *
import tree


def calculate_max_dot_product(state):
    state = state.as_single_qubit()
    a, b = state.initial_secret_params

    initial_state = tree.SingleQubit([a, b], state.initial_secret_params, state.resource_and_basis_parameters)
    initial_state.subs_b_by_a()

    state.subs_resource_and_basis_params()
    state = state.apply_basic_unitary()
    # state.pprint()

    optimized_state = state.apply_advanced_unitary()
    # optimized_state.pprint()

    params = state.subs_secret_params()
    print params

    initial_state.subs_params(params, with_val=True)
    dp = initial_state.dot_product(state)
    debug(dp=dp)

    optimized_state.subs_b_by_a()
    optimized_state.subs_params(params, with_val=True)
    print initial_state.dot_product(optimized_state)


def calculate_max_dot_product_new(state):
    state = state.as_single_qubit()
    a, b = state.initial_secret_params

    # Initial secret a|00\ + b|11\
    initial_state = tree.SingleQubit([a, b], state.initial_secret_params, state.resource_and_basis_parameters)
    initial_state.subs_b_by_a()

    state.subs_resource_and_basis_params()
    # state = state.apply_basic_unitary()
    # state.pprint()

    optimized_state = state.apply_advanced_unitary()
    # optimized_state.pprint()

    params = state.subs_secret_params()
    print params

    initial_state.subs_params(params, with_val=True)
    dp = initial_state.dot_product(state)
    debug(dp=dp)

    optimized_state.subs_b_by_a()
    optimized_state.subs_params(params, with_val=True)
    print initial_state.dot_product(optimized_state)


def main():
    root = tree.create_root()
    states = root.get_descendants(3)
    # state = states[random.randint(0, len(states) - 1)]
    # state = states[8]
    for  state in states:
        state.pprint()
        # calculate_max_dot_product_new(state)
        # print "="*15
        calculate_max_dot_product(state)

if __name__ == "__main__":
    main()

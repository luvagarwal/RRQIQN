from __future__ import division
import random
from sympy import sqrt

from utils import *
import tree


def calculate_max_dot_product(state):
    state = state.as_single_qubit()
    a, b = state.initial_secret_params

    initial_state = tree.SingleQubit([a, b], state.initial_secret_params, state.resouce_and_basis_params)
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


def main():
    root = tree.create_root()
    states = root.get_descendants(3)
    state = states[random.randint(0, len(states) - 1)]
    state.pprint()
    return calculate_max_dot_product(state)

if __name__ == "__main__":
    main()

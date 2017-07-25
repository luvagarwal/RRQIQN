from __future__ import division
import random
from sympy import sqrt

from utils import *
import tree
from paramvalues import ParamsSingleton


def drange(start, stop, step):
    r = start
    eps = 0.000001
    while r < stop - eps:
        yield r
        r += step

def analyze():
    root = tree.create_root()
    states = root.get_descendants(2)

    for two_qubit_state in states:
        entries = []
        repr_two_qubit = repr(two_qubit_state)
        print repr_two_qubit
        for pval in drange(0.1, 1, 0.1):
            for lval in drange(0.1, 1, 0.1):
                for qval in drange(0.1, 1, 0.1):
                    for aval in drange(0.1, 1, 0.1):
                        state = two_qubit_state.as_single_qubit()
                        a, b = state.initial_secret_params
                        p, l, q = state.resource_and_basis_parameters
                        ParamsSingleton().set_params_values({
                            p: pval,
                            l: lval,
                            q: qval,
                            a: aval, 
                        })
                        x, a_dash, dp_basic, dp_direct, dp_two_step = process_single_iteration(state)
                        entry = (format(pval), format(lval), format(qval), format(aval), format(x, upto=3), format(a_dash, upto=3), format(dp_basic), format(dp_direct))
                        # entry = ("x", x, "a_dash", a_dash, dp_basic, dp_direct)
                        entries.append(entry)
                print len(entries)

        with open('output', 'a') as f:
            entries = map(lambda entry: ", ".join(entry), entries)
            output = repr_two_qubit + "\n"
            output += "\n".join(entries) 
            output += "\n"
            f.write(output)


def process_single_iteration(state):
    a, b = state.initial_secret_params

    initial_state = tree.SingleQubit([a, b], state.initial_secret_params, state.resource_and_basis_parameters)
    initial_state.subs_b_by_a()

    state.subs_resource_and_basis_params()
    unitary_direct, optimized_state_direct = state.apply_advanced_unitary(return_unitary=True)

    received_network_state = state
    state = state.apply_basic_unitary()
    optimized_state = state.apply_advanced_unitary()
    # optimized_state.pprint()

    params = state.subs_secret_params()
    received_network_state.subs_params(params, with_val=True)
    a_dash = received_network_state.coeffs[0]

    initial_state.subs_params(params, with_val=True)
    # print initial_state.coeffs
    # debug(dp=dp)

    state.subs_b_by_a()
    state.subs_params(params, with_val=True)
    # print state.coeffs[0]

    optimized_state.subs_b_by_a()
    optimized_state.subs_params(params, with_val=True)
    # debug(optimized_dp=initial_state.dot_product(optimized_state))

    optimized_state_direct.subs_b_by_a()
    optimized_state_direct.subs_params(params, with_val=True)
    # print optimized_state_direct.coeffs
    # debug(optimized_dp_direct=initial_state.dot_product(optimized_state_direct))
    return unitary_direct, a_dash, initial_state.dot_product(state), initial_state.dot_product(optimized_state_direct), initial_state.dot_product(optimized_state)


def calculate_max_dot_product(state):
    state = state.as_single_qubit()
    a, b = state.initial_secret_params

    initial_state = tree.SingleQubit([a, b], state.initial_secret_params, state.resource_and_basis_parameters)
    initial_state.subs_b_by_a()

    state.subs_resource_and_basis_params()
    optimized_state_direct = state.apply_advanced_unitary()
    state = state.apply_basic_unitary()
    optimized_state = state.apply_advanced_unitary()
    # optimized_state.pprint()

    params = state.subs_secret_params()
    print params

    initial_state.subs_params(params, with_val=True)
    dp = initial_state.dot_product(state)
    debug(dp=dp)

    optimized_state.subs_b_by_a()
    optimized_state.subs_params(params, with_val=True)
    debug(optimized_dp=initial_state.dot_product(optimized_state))

    optimized_state_direct.subs_b_by_a()
    optimized_state_direct.subs_params(params, with_val=True)
    debug(optimized_dp_direct=initial_state.dot_product(optimized_state_direct))


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
    # main()
    analyze()



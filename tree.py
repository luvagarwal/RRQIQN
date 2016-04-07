#!/usr/bin/env python
from __future__ import division
from sympy import *
from sympy.parsing.sympy_parser import parse_expr

class node:
    def __init__(self, level, typ, coeff, coeff_comm=1, max_level=4):
        self.typ = typ
        self.coeff = coeff
        self.coeff_comm = coeff_comm
        self.level = level
        self.max_level = max_level
        if(self.level < self.max_level):
            self.bacche = self.find_bacche()
    
    def find_bacche(self):
        a, b = self.coeff
        c = self.coeff_comm
        level = self.level
        # p, l = symbols('p%s l%s'%(level+1, level+1))
        p, l, q, Q, L, P = symbols('p l q Q L P')
        if self.typ == 1:
            return  [        
                node(level+1, 1, [a, b*q*l.conjugate()], c*Q*L),
                node(level+1, 1, [a*l, -b*q], c*Q*L),
                node(level+1, 2, [a*q, b*p.conjugate()], c*Q*P),
                node(level+1, 2, [a*q*p, -b], c*Q*P)
            ]
        else: # type = 2 -> a|01\ + b|10\
            return  [
                node(level+1, 2, [b*q*l.conjugate(), a], c*Q*L),
                node(level+1, 2, [-b*q, a*l], c*Q*L),
                node(level+1, 1, [b*p.conjugate(), a*q], c*Q*P),
                node(level+1, 1, [-b, a*q*p], c*Q*P)
            ]


def compact(a):
    return collect(a, a.atoms(Symbol), func=factor)

def parse_level(node, level, t="phi_plus"):
    d = {
        'phi_plus': 0,
        'phi_minus': 1, 
        'psi_plus': 2,
        'psi_minus': 3,
        }
    if node.level == level-1:
        ans = [0]*4
        n = node.bacche[d[t]]
        if n.typ == 1:
            ans[0] += n.coeff[0]
            ans[1] += n.coeff[1]
        if n.typ == 2:
            ans[2] += n.coeff[0]
            ans[3] += n.coeff[1]
        return ans

    ans = [0]*4
    for n in node.bacche:
        a = parse_level(n, level, t)
        ans = [x + y for x, y in zip(a, ans)]
    return ans

def get_node():
    a, b = symbols('a b')
    n = node(0, 1, [a, b])
    return n 

def print_level(level, t="phi_plus"):
    n = get_node()
    ans = parse_level(n, level, t)
    return ans

def _nodes_as_array(n, level):
    " return the value of nodes at different level in array form "
    if n.level == level - 1:
        return n.bacche
    ans = []
    for node in n.bacche:
        ans += _nodes_as_array(node, level)
    return ans

def nodes_as_array(level, typ=None):
    n = get_node()
    ans = _nodes_as_array(n, level)
    f_ans = []
    for a in ans:
        if a.typ == 1:
            if typ and typ != a.typ:
                continue
            tmp = "%s*(%s*Symbol('|00|') + %s*Symbol('|11|'))"%(a.coeff_comm, a.coeff[0], a.coeff[1])
            tmp = parse_expr(tmp)
            f_ans.append(tmp)
        else:
            if typ and typ != a.typ:
                continue
            tmp = "%s*(%s*Symbol('|01|') + %s*Symbol('|10|')"%(a.coeff_comm, a.coeff[0], a.coeff[1])
            tmp = parse_expr(tmp)
            f_ans.append(tmp)
    return f_ans

def main():
    pass

if __name__ == "__main__":
    nodes_as_array(3)

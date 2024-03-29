r"""
Homogeneous differential polynomial equation solver
"""
from collections import defaultdict
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.misc.verbose import verbose
from .differential_polynomial_ring import DifferentialPolynomial

# TODO: system of equations, e.g. for H

def solve_homogeneous_diffpoly(target, source, unknowns):
    """
    Return a solution of a homogeneous differential polynomial equation.

    INPUT:

    - ``target`` -- a homogeneous differential polynomial, the right-hand side of the equation

    - ``source`` -- a homogeneous differential polynomial, the left-hand side of the equation

    - ``unknowns`` -- a list of fibre variables, such that the total derivatives of those variables appear in ``source``

    ALGORITHM:

    Builds an ansatz based on the homogeneity, and solves the arising linear system.
    """
    if not isinstance(target, DifferentialPolynomial):
        raise ValueError('target must be a differential polynomial')
    if not isinstance(source, DifferentialPolynomial):
        raise ValueError('source must be a differential polynomial')

    target_monomials = set(target.monomials())
    if not target.is_fibre_degree_homogeneous():
        raise ValueError('target must be fibre degree homogeneous')
    target_degrees = next(iter(target_monomials)).fibre_degrees()
    verbose('target degrees: {}'.format(target_degrees), level=1)

    if not target.is_weight_homogeneous():
        raise ValueError('target must be weight homogeneous')
    target_weights = next(iter(target_monomials)).weights()
    verbose('target weights: {}'.format(target_weights), level=1)

    R = target.parent()
    unknowns_derivatives = defaultdict(set)
    all_unknowns = set(unknowns)
    for v in source.variables():
        if v in R.jet_variables():
            w, s = v.variable_subscript()
            if w in unknowns:
                unknowns_derivatives[w].add(v)
                all_unknowns.add(v)

    source_part = {}
    for v in unknowns:
        subs = {v : R.zero() for v in all_unknowns}
        subs[v] = v
        for w in unknowns_derivatives[v]:
            subs[w] = w
        source_part[v] = source.subs(subs)

    # infer weights and degrees of ansatz using source and target
    ansatz_degrees = defaultdict(set)
    ansatz_weights = defaultdict(set)
    for v in all_unknowns:
        w, s = v.variable_subscript()
        for m in (source_part[w] // v).monomials():
            ansatz_degree = target_degrees - m.fibre_degrees() # equation is linear in unknowns
            ansatz_degree.set_immutable()
            ansatz_degrees[w].add(ansatz_degree)
            ansatz_weight = target_weights - m.weights() - v.weights()
            ansatz_weight.set_immutable()
            ansatz_weights[w].add(ansatz_weight)

    verbose('ansatz degrees: {}'.format(dict(ansatz_degrees)), level=1)
    verbose('ansatz weights: {}'.format(dict(ansatz_weights)), level=1)

    # generate monomials. TODO: optional?
    ansatz_monomials = defaultdict(set)
    for v in unknowns:
        for d in ansatz_degrees[v]:
            for w in ansatz_weights[v]:
                ansatz_monomials[v].update(R.homogeneous_monomials(d, w, max_differential_orders=R.max_differential_orders()))

    verbose('ansatz #monomials: ' + str({d : len(m) for d,m in ansatz_monomials.items()}), level=1)
    verbose('ansatz monomials: {}'.format(ansatz_monomials), level=2)

    # update target monomials by substituting ansatz for fibre variables (and total derivatives)
    for v in unknowns:
        pre_subs = dict(zip(all_unknowns, [R.zero()]*len(all_unknowns)))
        for m in ansatz_monomials[v]:
            subs = pre_subs.copy()
            subs[v] = m
            admissible = True
            for w in unknowns_derivatives[v]:
                _, s = w.variable_subscript()
                try:
                    subs[w] = m.total_derivative(*s)
                except ValueError:
                    admissible = False
                    break
            if not admissible:
                continue
            f = source_part[v].subs(subs)
            target_monomials.update(f.monomials())

    target_basis = list(target_monomials)

    verbose('len(target_basis) == {}'.format(len(target_basis)), level=1)
    verbose('target basis: {}'.format(target_basis), level=2)

    M = matrix(R.base_ring(), len(target_basis), 0, sparse=True)
    ansatz_basis = []
    for v in unknowns:
        pre_subs = dict(zip(all_unknowns, [R.zero()]*len(all_unknowns)))
        for m in ansatz_monomials[v]:
            subs = pre_subs.copy()
            subs[v] = m
            admissible = True
            for w in unknowns_derivatives[v]:
                _, s = w.variable_subscript()
                try:
                    subs[w] = m.total_derivative(*s)
                except ValueError:
                    admissible = False
                    break
            if not admissible:
                continue
            ansatz_basis.append((v, m))
            f = source_part[v].subs(subs)
            V = vector(R.base_ring(), len(target_basis), {target_basis.index(m) : c for c,m in f}, sparse=True)
            M = M.augment(V)

    verbose('len(ansatz_basis) == {}'.format(len(ansatz_basis)), level=1)
    verbose('ansatz basis: {}'.format(ansatz_basis), level=2)

    target_vector = vector(R.base_ring(), len(target_basis), {target_basis.index(m) : c for c,m in target}, sparse=True)

    solution_vector = M.solve_right(target_vector)
    solution = {v : R.zero() for v in unknowns}
    for i in range(len(ansatz_basis)):
        v, b = ansatz_basis[i]
        solution[v] += b*solution_vector[i]
    return solution

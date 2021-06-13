from .differential_polynomial_ring import DifferentialPolynomial
from collections import defaultdict
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.misc.verbose import verbose

# TODO: system of equations, e.g. for H

def solve_homogeneous_diffpoly(target, source, unknowns):
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

    verbose('ansatz monomials: {}'.format(ansatz_monomials), level=1)

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
                    subs[w] = m.tdiff(*s)
                except ValueError:
                    admissible = False
                    break
            if not admissible:
                continue
            f = source_part[v].subs(subs)
            target_monomials.update(f.monomials())

    target_basis = list(target_monomials)

    verbose('target basis: {}'.format(target_basis), level=1)

    M = matrix(R.base_ring(), len(target_basis), 0)
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
                    subs[w] = m.tdiff(*s)
                except ValueError:
                    admissible = False
                    break
            if not admissible:
                continue
            ansatz_basis.append((v, m))
            f = source_part[v].subs(subs)
            V = vector(R.base_ring(), [f.monomial_coefficient(b) for b in target_basis])
            M = M.augment(V)

    verbose('ansatz basis: {}'.format(ansatz_basis), level=1)

    target_vector = vector(R.base_ring(), [target.monomial_coefficient(b) for b in target_basis])

    solution_vector = M.solve_right(target_vector)
    solution = {v : R.zero() for v in unknowns}
    for i in range(len(ansatz_basis)):
        v, b = ansatz_basis[i]
        solution[v] += b*solution_vector[i]
    return solution

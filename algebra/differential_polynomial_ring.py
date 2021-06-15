from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.modules.free_module_element import vector
from sage.rings.integer_ring import ZZ
from sage.combinat.integer_vector import IntegerVectors
from sage.misc.misc_c import prod
from sage.symbolic.expression import is_Expression
from sage.calculus.var import var, function
from itertools import product, combinations
from util.jet_variables import SubstituteJetVariables, SubstituteTotalDerivatives

class DifferentialPolynomial:
    def __init__(self, parent, polynomial):
        self._parent = parent
        if not polynomial in self._parent._polynomial_ring:
            raise ValueError('polynomial must be in polynomial ring of the parent')
        self._polynomial = polynomial
        
    def __repr__(self):
        return repr(self._polynomial)
    
    def parent(self):
        return self._parent
    
    def __eq__(self, other):
        if other == 0:
            return self._polynomial.is_zero()
        return self.parent() == other.parent() and self._polynomial == other._polynomial
    
    def __ne__(self, other):
        return not self.__eq__(other)
    
    def __hash__(self):
        return hash((self._parent, self._polynomial))
    
    def copy(self):
        return __class__(self._parent, self._polynomial)
    
    __pos__ = copy
    
    def __neg__(self):
        return __class__(self._parent, -self._polynomial)
    
    def __add__(self, other):
        if isinstance(other, __class__):
            return __class__(self._parent, self._polynomial + other._polynomial)
        else:
            return __class__(self._parent, self._polynomial + other)
    
    def __radd__(self, other):
        return self + other
    
    def __sub__(self, other):
        if isinstance(other, __class__):
            return __class__(self._parent, self._polynomial - other._polynomial)
        else:
            return __class__(self._parent, self._polynomial - other)
    
    def __rsub__(self, other):
        return -(self - other)
    
    def __mul__(self, other):
        if isinstance(other, __class__):
            return __class__(self._parent, self._polynomial * other._polynomial)
        else:
            return __class__(self._parent, self._polynomial * other)
    
    def __rmul__(self, other):
        return self * other
    
    def __pow__(self, other):
        return __class__(self._parent, self._polynomial^other)
    
    def __truediv__(self, other):
        return __class__(self._parent, self._polynomial / other)

    def __floordiv__(self, other):
        if isinstance(other, __class__):
            return __class__(self._parent, self._polynomial // other._polynomial)
        else:
            return __class__(self._parent, self._polynomial // other)
    
    def __mod__(self, other):
        if isinstance(other, __class__):
            return __class__(self._parent, self._polynomial % other._polynomial)
        else:
            return __class__(self._parent, self._polynomial % other)
    
    def is_zero(self):
        return self._polynomial.is_zero()

    def _derivative_once(self, x):
        return __class__(self._parent, self._polynomial.derivative(x._polynomial))
    
    def derivative(self, *x):
        result = self.copy()
        for v in x:
            result = result._derivative_once(v)
        return result
    
    diff = derivative
    
    def _total_derivative_once(self, x):
        result = self._polynomial.derivative(x._polynomial)
        for v in self._polynomial.variables():
            # TODO: check fibre variable
            result += self._polynomial.derivative(v) * self._parent._diff_single_var(v, x._polynomial)
        return __class__(self._parent, result)
    
    def total_derivative(self, *x):
        result = self.copy()
        for v in x:
            result = result._total_derivative_once(v)
        return result
    
    tdiff = total_derivative
    
    def _integral_monomials_once(self, x):
        # TODO: assumes monomial
        result = set([])
        for v in self._polynomial.variables():
            f = (self._polynomial // v) * self._parent._integrate_single_var(v, x._polynomial)
            result.add(__class__(self._parent, f))
        return result
    
    def integral_monomials(self, *x):
        result = set(self.monomials())
        for v in x:
            # TODO: improve?
            result2 = set([])
            for m in result:
                result2.update(m._integral_monomials_once(v))
            result = result2
        return result
    
    def substitute(self, arg):
        if not isinstance(arg, dict):
            raise ValueError('can only substitute dict')
        poly_arg = {k._polynomial : v._polynomial for (k,v) in arg.items()}
        return __class__(self._parent, self._polynomial.subs(poly_arg))
    
    subs = substitute

    def _symbolic_(self, ring):
        result = ring(self._polynomial)
        result = self._parent._subs_tot_ders(result)
        return result
    
    def degree(self):
        return self._polynomial.degree()
    
    def variables(self):
        return tuple(__class__(self._parent, v) for v in self._polynomial.variables())
    
    def monomials(self):
        return tuple(__class__(self._parent, m) for m in self._polynomial.monomials())
    
    def coefficients(self):
        return self._polynomial.coefficients()
    
    def monomial_coefficient(self, m):
        return self._polynomial.monomial_coefficient(m._polynomial)
    
    def __iter__(self):
        for (c,m) in self._polynomial:
            yield (c, __class__(self._parent, m))
    
    def variable_subscript(self):
        idx = self._parent._var_to_idx[self._polynomial]
        fibre_idx = idx[0]
        subscript_idx = idx[1:]
        base_vars = self._parent.base_variables()
        return self._parent.fibre_variable(fibre_idx), tuple(sum(([base_vars[i]]*subscript_idx[i] for i in range(len(subscript_idx))), []))
    
    def weights(self):
        base_dim = self._parent.base_dim()
        w = [0 for i in range(base_dim)]
        mon = self.monomials()
        if len(mon) != 1:
            raise ValueError('weights are only defined for monomials')
        u = self._parent.gens()
        e = mon[0]._polynomial.exponents()[0]
        for i in range(base_dim, len(e)):
            w_i = self._parent._single_var_weights(u[i]._polynomial)
            w_i = [w_ij*e[i] for w_ij in w_i]
            w = [w_a + w_b for (w_a,w_b) in zip(w,w_i)]
        return vector(ZZ, w, immutable=True)
    
    def is_weight_homogeneous(self):
        return len(set(m.weights() for m in self.monomials())) == 1

    def fibre_degrees(self):
        mon = self.monomials()
        if len(mon) != 1:
            raise ValueError('fibre degrees are only defined for monomials')
        e = mon[0]._polynomial.exponents()[0]
        base_dim = self._parent.base_dim()
        fibre_vars = self._parent.fibre_variables()
        fibre_dim = len(fibre_vars)
        degrees = [0 for i in range(fibre_dim)]
        for i in range(base_dim, len(e)):
            if e[i] == 0:
                continue
            v = self._parent.gen(i)
            if i < base_dim + fibre_dim:
                w = v
            else:
                w, _ = v.variable_subscript()
            degrees[fibre_vars.index(w)] += e[i]
        return vector(ZZ, degrees, immutable=True)
    
    def is_fibre_degree_homogeneous(self):
        return len(set(m.fibre_degrees() for m in self.monomials())) == 1

class DifferentialPolynomialRing:
    element_class = DifferentialPolynomial
    
    def __init__(self, base_ring, fibre_names, base_names, max_differential_orders):
        self._fibre_names = tuple(fibre_names)
        self._base_names = tuple(base_names)
        self._max_differential_orders = tuple(max_differential_orders)
        base_dim = len(self._base_names)
        fibre_dim = len(self._fibre_names)
        jet_names = []
        idx_to_name = {}
        for fibre_idx in range(fibre_dim):
            u = self._fibre_names[fibre_idx]
            idx_to_name[(fibre_idx,) + tuple([0]*base_dim)] = u
            for d in range(1, max_differential_orders[fibre_idx]+1):
                for multi_index in IntegerVectors(d, base_dim):
                    v = '{}_{}'.format(u, ''.join(self._base_names[i]*multi_index[i] for i in range(base_dim)))
                    jet_names.append(v)
                    idx_to_name[(fibre_idx,) + tuple(multi_index)] = v
        self._polynomial_ring = PolynomialRing(base_ring, base_names + fibre_names + tuple(jet_names))
        self._idx_to_var = {idx : self._polynomial_ring(idx_to_name[idx]) for idx in idx_to_name}
        self._var_to_idx = {jet : idx for (idx,jet) in self._idx_to_var.items()}
        # for conversion:
        base_vars = [var(b) for b in self._base_names]
        symbolic_functions = [function(f)(*base_vars) for f in self._fibre_names]
        self._subs_jet_vars = SubstituteJetVariables(symbolic_functions)
        self._subs_tot_ders = SubstituteTotalDerivatives(symbolic_functions)
    
    def __repr__(self):
        return 'Differential Polynomial Ring in {} over {}'.format(', '.join(map(repr, self._polynomial_ring.gens())), self._polynomial_ring.base_ring())

    def base_ring(self):
        return self._polynomial_ring.base_ring()
    
    def _first_ngens(self, n):
        return tuple(self.element_class(self, self._polynomial_ring.gen(i)) for i in range(n))
    
    def gens(self):
        return self._first_ngens(self._polynomial_ring.ngens())
    
    def gen(self, i):
        return self.element_class(self, self._polynomial_ring.gen(i))

    def base_variables(self):
        return self._first_ngens(len(self._base_names))
    
    def base_dim(self):
        return len(self._base_names)
    
    def fibre_variable(self, i):
        return self.element_class(self, self._polynomial_ring.gen(len(self._base_names) + i))
    
    def fibre_variables(self):
        base_dim = len(self._base_names)
        fibre_dim = len(self._fibre_names)
        return tuple(self.element_class(self, self._polynomial_ring.gen(base_dim + i)) for i in range(fibre_dim))
    
    def fibre_dim(self):
        return len(self._fibre_names)
    
    def jet_variables(self):
        base_dim = len(self._base_names)
        fibre_dim = len(self._fibre_names)
        whole_dim = self._polynomial_ring.ngens()
        return tuple(self.element_class(self, self._polynomial_ring.gen(i)) for i in range(base_dim + fibre_dim, whole_dim))
    
    def max_differential_orders(self):
        return self._max_differential_orders
    
    def _single_var_weights(self, u):
        return self._var_to_idx[u][1:]
    
    def _diff_single_var(self, u, x):
        x_idx = self._polynomial_ring.gens().index(x)
        u_idx = self._var_to_idx[u]
        du_idx = list(u_idx)
        du_idx[1 + x_idx] += 1
        du_idx = tuple(du_idx)
        if du_idx in self._idx_to_var:
            return self._idx_to_var[du_idx]
        else:
            raise ValueError("can't differentiate {} any further with respect to {}".format(u, x))
    
    def _integrate_single_var(self, u, x):
        x_idx = self._polynomial_ring.gens().index(x)
        u_idx = self._var_to_idx[u]
        if u_idx[1 + x_idx] == 0:
            raise ValueError("can't integrate {} any further with respect to {}".format(u,x))
        iu_idx = list(u_idx)
        iu_idx[1 + x_idx] -= 1
        iu_idx = tuple(iu_idx)
        return self._idx_to_var[iu_idx]
            
    def __call__(self, arg):
        if is_Expression(arg):
            arg = self._subs_jet_vars(arg)
        return self.element_class(self, self._polynomial_ring(arg))
    
    def zero(self):
        return self.element_class(self, self._polynomial_ring.zero())
    
    def one(self):
        return self.element_class(self, self._polynomial_ring.one())

    def homogeneous_monomials(self, fibre_degrees, weights, max_differential_orders=None):
        fibre_vars = self.fibre_variables()
        if not len(fibre_degrees) == len(fibre_vars):
            raise ValueError('length of fibre_degrees vector must match number of fibre variables')
        base_vars = self.base_variables()
        if not len(weights) == len(base_vars):
            raise ValueError('length of weights vector must match number of base variables')
        monomials = []
        fibre_degree = sum(fibre_degrees)
        fibre_indexes = {}
        fibre_idx = 0
        for i in range(len(fibre_degrees)):
            for j in range(fibre_degrees[i]):
                fibre_indexes[fibre_idx] = i
                fibre_idx += 1
        proto = sum([[fibre_vars[i]]*fibre_degrees[i] for i in range(len(fibre_degrees))], [])
        for V in product(*[IntegerVectors(w, fibre_degree) for w in weights]):
            total_differential_order = [0 for i in range(fibre_degree)]
            term = [p for p in proto]
            skip = False
            for j in range(fibre_degree):
                fibre_idx = fibre_indexes[j]
                for i in range(len(base_vars)):
                    if V[i][j] > 0:
                        total_differential_order[j] += V[i][j]
                        if max_differential_orders is not None and total_differential_order[j] > max_differential_orders[fibre_idx]:
                            skip = True
                            break
                        term[j] = term[j].tdiff(*([base_vars[i]]*V[i][j]))
                if skip:
                    break
            if not skip:
                monomials.append(prod(term))
        return monomials
    
def TD(a,*x):
    return a.tdiff(*x)

from sage.symbolic.expression_conversions import ExpressionTreeWalker
from sage.symbolic.ring import SR

class SubstituteJetVariables(ExpressionTreeWalker):
    def __init__(self, functions):
        self._functions = [f.operator() for f in functions]
        self._base_vars = [f.operands() for f in functions]

    def composition(self, ex, operator):
        if operator in self._functions:
            return SR.var(str(operator))
        else:
            return super().composition(ex, operator)

    def derivative(self, ex, operator):
        f = operator.function()
        if f in self._functions:
            f_idx = self._functions.index(f)
            base_vars = self._base_vars[f_idx]
            return SR.var(str(operator.function()) + '_' + ''.join(str(base_vars[i]) for i in operator.parameter_set()))
        else:
            return operator(*[self(_) for _ in ex.operands()])

class SubstituteTotalDerivatives(ExpressionTreeWalker):
    def __init__(self, functions):
        self._functions = functions
        self._function_names = [str(f.operator()) for f in functions]
        self._base_vars = [f.operands() for f in functions]
        self._base_vars_names = [[str(v) for v in x] for x in self._base_vars]

    def symbol(self, x):
        str_x = str(x)
        if str_x in self._function_names:
            f_idx = self._function_names.index(str_x)
            return self._functions[f_idx]
        elif '_' in str_x:
            parts = str_x.split('_')
            if len(parts) == 2:
                variable, subscripts = parts
                if variable in self._function_names:
                    f_idx = self._function_names.index(variable)
                    base_vars_names = self._base_vars_names[f_idx]
                    base_vars = self._base_vars[f_idx]
                    result = self._functions[f_idx]
                    for v in subscripts:
                        v_idx = base_vars_names.index(v)
                        result = result.derivative(base_vars[v_idx])
                    return result
                else:
                    return super().symbol(x)
            else:
                return super().symbol(x)
        else:
            return super().symbol(x)

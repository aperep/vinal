import sympy

class Cone:

    def initial_pair(self):
        pivot_rows = self.A_matrix().pivot_rows()
        A0 = [self.A()[pivot] for pivot in pivot_rows]
        Ac = [self.A()[i] for i in range(len(self.A())) if i not in pivot_rows]
        from sage.matrix.constructor import identity_matrix, matrix
        I = identity_matrix(self.base_ring(), self.dim())
        R = matrix(self.base_ring(), A0).solve_right(I)
        return self.pair_class(self, A0, R.columns()), list(Ac)
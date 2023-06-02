import dinv
import sympy as sp

if __name__ == '__main__':
    example = sp.Matrix([
        [1, 0, 0],
        [-1, 1, 0],
        [0, -1, 0],
        [0, 0, 1],
        [0, 0, -1],
    ])

    print(f'd-invariant table for lattice:\n{sp.pretty(example)}\n')
    table = dinv.dinv_table(example)
    print(table)

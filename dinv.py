import itertools
import sympy as sp
import numpy as np
from tabulate import tabulate

# Compute the d-invariant of a lattice, assuming it is the lattice of integer flows or cuts of
# a graph.


def dinv_table(basis):
    n, m = basis.shape

    # Compute the matrix that projects Z^n onto the lattice.
    project = sp.zeros(n)

    for j in range(m):
        vec = basis.col(j)
        project += 1 / (vec.dot(vec)) * vec @ vec.transpose()

    # Compute a transformation from Z^n to the basis. The pseudoinverse is such a transformation.
    change = (basis.transpose() @ basis).inv() @ basis.transpose()

    transform = change @ project

    # Compute short covectors/orientation covectors.
    orientation_covectors = map(
        lambda coeffs: sp.matrices.Matrix(coeffs),
        itertools.product(
            (1, -1),
            repeat=n
        )
    )

    table = []
    for o in orientation_covectors:
        row = []

        chi = transform @ o
        chi_euclidean = project @ o
        normsq_chi = chi_euclidean.dot(chi_euclidean)
        d = (normsq_chi - m) / 4

        row.append(str(tuple(o.transpose())))
        row.append(str(tuple(chi.transpose())))
        row.append(str(normsq_chi))
        row.append(str(d))

        table.append(row)

    return tabulate(table, headers=['Short', 'Chi', '|Chi|', 'd'], numalign='left')

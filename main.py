import sympy as sp
import networkx as nx
import typing
import itertools

# Set up printing.
sp.init_printing()
old_print = print
pprint_type = sp.MatrixBase | typing.List | typing.Tuple | typing.Set |\
    typing.Dict


def my_print(obj):
    if isinstance(obj, pprint_type):
        sp.pprint(obj)
    else:
        old_print(obj)


print = my_print

# Take edge list as input.
edge_list = sp.matrices.Matrix(
    [
        [0, 1],
        [0, 1],
        [1, 0],
    ],
)

n_edges = edge_list.shape[0]
n_vertices = max(sp.utilities.iterables.flatten(edge_list)) + 1

# Compute the networkx MultiGraph from the edge list.
graph = nx.MultiDiGraph()
graph.add_nodes_from(range(n_vertices))
for i in range(n_edges):
    u, v = edge_list.row(i)
    graph.add_edge(u, v, key=i)


# Compute the flow_basis_matrix.
#
# Matrix indexed by (i, j) where the jth column is the column vector of a flow
# written in the edge-basis.
spanning_tree = nx.maximum_spanning_tree(graph.to_undirected(as_view=True))
spanning_edges = set()
for _u, _v, edge in spanning_tree.edges(keys=True):
    spanning_edges.add(edge)
not_tree_edges = set(range(n_edges)) - spanning_edges
flow_basis = []
for starting_edge in not_tree_edges:
    circuit = sp.matrices.zeros(n_edges, 1)
    circuit[starting_edge] = 1
    u, v = edge_list.row(starting_edge)
    assert graph.has_edge(u=u, v=v, key=starting_edge)
    spath = nx.shortest_path(spanning_tree, source=v, target=u)
    for u, v in itertools.pairwise(spath):
        edges_in_st = spanning_tree.get_edge_data(u, v)
        assert len(edges_in_st) == 1
        edge = list(edges_in_st)[0]
        if graph.has_edge(u, v, key=edge):
            circuit[edge] = 1
        elif graph.has_edge(v, u, key=edge):
            circuit[edge] = -1
        else:
            raise ValueError('Edge not found in graph.')
    flow_basis.append(circuit)
flow_basis_matrix = sp.matrices.Matrix([flow_basis])
print('flow lattice basis matrix:')
print(flow_basis_matrix)


def dual_basis_matrix(basis_matrix):
    pseudo_inverse = (basis_matrix.H * basis_matrix).inv() * basis_matrix.H
    return pseudo_inverse.H


fdual_basis_matrix = dual_basis_matrix(flow_basis_matrix)
print('dual of flow lattice basis matrix:')
print(fdual_basis_matrix)


# Calculate the short vectors of the lattice of integer flows via the
# surjection from the short vectors of the unimdoular lattice.
short_unimod = map(
    lambda coeffs: sp.matrices.Matrix(coeffs),
    itertools.product([1, -1], repeat=n_edges)
)
short_flow = map(
    lambda v: fdual_basis_matrix.H * v,
    short_unimod
)
distinct_short_flow = set()
for v in short_flow:
    v = sp.matrices.ImmutableMatrix(v)
    # Check v not in same equiv. class as any other vector u so far by checking
    # coefficients of (u - v).
    if not any(all(c % 2 == 0 for c in (v - u)) for u in distinct_short_flow):
        distinct_short_flow.add(v)
print("short vectors (mod 2 " + "\U0001D509" + "(G)):")
sp.pprint(set(distinct_short_flow))


# Calculate the d-invariant of the graph.
flow_rank = flow_basis_matrix.rank()
# In the flow basis, the new inner product matrix is given by:
flow_gram_matrix = flow_basis_matrix.H * flow_basis_matrix
d_map = {
    x:
    (sp.MatrixBase.dot(x.H * flow_gram_matrix, x) - flow_rank) / 4
    for x in distinct_short_flow
}

print(d_map)

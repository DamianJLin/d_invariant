import numpy as np
import itertools

# Take edge list as input.
edge_list = np.asarray(
    [
        [0, 1],
        [0, 1],
        [0, 1],
    ],
    dtype=int
)
edge_list = np.asarray(
    [
        [0, 3],
        [2, 0],
        [0, 1],
        [2, 1],
        [1, 3],
        [2, 3],
    ],
    dtype=int
)
edge_list = np.asarray(
    [
        [0, 1],
        [1, 2],
        [0, 2],
    ],
    dtype=int
)

# Compute the adjacency matrix.
#
# Matrix indexed (i, j) whose (i, j) the entry represents the jth edge's
# coefficient of the boundary map in the ith vertex.
n_edges = edge_list.shape[0]
n_vertices = np.max(edge_list.flatten() + 1)
adjacency_matrix = np.zeros((n_vertices, n_edges), dtype=int)
for j, (u, v) in enumerate(edge_list):
    adjacency_matrix[u, j] = -1
    adjacency_matrix[v, j] = 1


# Compute circuits of the adjacency matrix, corresponding to circuits of the
# graph.
#
# TODO: Implelemt efficient algorithm (currently factorial time, from a brief
# survey, there are known poly time algorithms).
edges = frozenset(range(n_edges))
edges_powerset = set(
    itertools.chain.from_iterable(
            itertools.combinations(edges, r) for r in range(len(edges) + 1)
    )
)
edges_powerset = {frozenset(es) for es in edges_powerset}
independents = {
    es for es in edges_powerset
    if len(es) == 0
    or np.linalg.matrix_rank(adjacency_matrix[np.ix_((0, 1), tuple(es))])
    == len(es)
}
dependents = edges_powerset.difference(independents)
circuits = {
    es for es in dependents
    if all(
        set(es) - {e}
        in independents
        for e in es)
}
# Sort two layers deep (sort each circuit, and sort circuits).
circuits = sorted(tuple(sorted(circuit)) for circuit in circuits)
n_circuits = len(circuits)
print(f'{circuits = }')


def boundary(edge_as_list):
    u, v = edge_as_list
    boundary = np.zeros((n_edges, 1), dtype=int)
    boundary[u, 0] = -1
    boundary[v, 0] = 1
    return boundary


# Compute the flow matrix.
#
# Matrix indexed by (i, j) where the jth column is the column vector of the jth
# circuit written in the edge-basis, with entries signs altered to make the
# circuit into a flow (ensure it is in the kernel of the boundary map).
flow_matrix = np.zeros((n_edges, n_circuits), dtype=int)
for j, circuit in enumerate(circuits):
    # k over edges by order in circuit, j over edges by order in edge_list.
    # Last edge, but weighted by the assigned coefficient.
    last_wedge = None
    for k, i in enumerate(circuit):
        edge = boundary(edge_list[i])
        if k == 0:
            w = 1
        else:
            w = -np.sign(np.dot(last_wedge.flatten(), edge.flatten()))
        flow_matrix[i, j] = w
        last_wedge = w * edge
        assert w == 1 or w == -1

# Check that the circuit matrix's columns are indeed in the kernel of the
# boundary map. TODO: Test and remove.
for j, col in enumerate(flow_matrix.transpose()):
    circuit = circuits[j]
    assert np.array_equal(
        sum(boundary(edge_list[i]) * col[i] for i in circuit),
        np.zeros((n_edges, 1), dtype=int)
    )
print(f'flow_matrix = \n{flow_matrix}')

# Compute the flow_basis_matrix.
#
# Matrix indexed by (i, j) where the jth column is the column vector of a flow
# written in the edge-basis. Since entries of flow_matrix necessarily have unit
# norm, they span the lattice of integer flows. Thus we only need to ensure
# independence.
r = np.linalg.matrix_rank(flow_matrix)
flow_basis_matrix = flow_matrix[:, :r]
print(f'flow_basis_matrix = \n{flow_basis_matrix}')

fdual_basis_matrix = np.linalg.pinv(flow_basis_matrix)
print(f'fdual_basis_matrix = \n{fdual_basis_matrix}')

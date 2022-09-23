import sys
import csv
import pathlib
import sympy as sp
import lib.dinv as dinv


if __name__ == '__main__':
    sp.init_printing()
    verbose = False

    # Parse args.
    args = sys.argv[1:]
    args_no_flags = []
    for arg in args:
        if arg[0] == '-':
            if arg == '-v':
                verbose = True
        else:
            args_no_flags.append(arg)
    args = args_no_flags

    if not args:
        raise ValueError("Argument(s) missing: edge list .csv file.")

    # Set up filepath.
    parent_path = pathlib.Path(__file__).parent
    graph_paths = [parent_path / fname for fname in args]

    # Compute d-invariant(s).
    for i, graph_path in enumerate(graph_paths):
        if not graph_path.is_file():
            raise FileNotFoundError(f"{graph_path} not found.")

        # Take edge list as input from graph.csv.
        with open(graph_path) as csvfile:
            reader = csv.reader(csvfile)
            edge_list = []
            for row in reader:
                if len(row) != 2:
                    raise ValueError("graph.csv must have row length 2.")
                if not all(entry.strip().isdigit() for entry in row):
                    raise ValueError("graph.csv entries must be integers.")
                edge_list.append([int(entry.strip()) for entry in row])
            if not edge_list:
                raise ValueError("graph.csv must not be an empty file.")

            edge_list = sp.matrices.Matrix(edge_list)
            chi, d = dinv.compute_d_invariant(edge_list)

            # TODO: Write to file.

            # Printing.
            if verbose:
                sp.pprint(f'd-invariant for graph in {graph_path.name}')
                sp.pprint("short vectors (mod 2 " + "\U0001D509" + "(G)):")
                sp.pprint(chi)
                sp.pprint('d mapping:')
                sp.pprint(d)
                if not i == len(args) - 1:
                    print('')

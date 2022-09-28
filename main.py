import sys
import pathlib
import sympy as sp
import pickle
import src.dinv as dinv
import src.tait as tait
import progress.bar as prog

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

    if not len(args) == 1:
        raise ValueError("Invalid arguments.")
    input_relpath = args[0]

    # Set up filepath.
    parent_path = pathlib.Path(__file__).parent
    input_path = parent_path / input_relpath
    output_path = parent_path / 'out'

    # Clear files in output_path.
    for f in output_path.glob('*'):
        f.unlink()

    # Set up progress bar.
    n_input = sum(1 for line in open(input_path))
    bar = prog.IncrementalBar(
        'Computing d-invariants:',
        max=n_input,
        color='green',
    )

    with open(input_path, 'r') as file_in:
        for line in file_in:
            data = line.rstrip().split(' ')
            name, gauss_code = data
            tnorm, tdual = tait.tait_graph_edge_lists(gauss_code)
            tnorm_dinv = dinv.compute_d_invariant(tnorm)
            tdual_dinv = dinv.compute_d_invariant(tdual)

            tait_graphs_path = output_path / f'{name}.tait'
            with open(tait_graphs_path, 'w') as file_out:
                for edge in tnorm:
                    file_out.write(str(edge) + '\n')
                file_out.write('\n')
                for edge in tdual:
                    file_out.write(str(edge) + '\n')

            dinv_path = output_path / f'{name}.dinv'
            with open(dinv_path, 'w') as file_out:
                file_out.write(sp.pretty((tnorm_dinv[1])))
                file_out.write('\n\n')
                file_out.write(sp.pretty((tdual_dinv[1])))

            dinv_pickle_path = output_path / f'{name}.dinv.pickle'
            with open(dinv_pickle_path, 'w') as file_out:
                pickle.dump(
                    (tnorm_dinv, tdual_dinv),
                    file_out,
                    protocol=pickle.HIGHEST_PROTOCOL
                )

            bar.next()
        bar.finish()

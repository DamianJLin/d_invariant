# dinv
Script written to compute the d-invariant of a knot/virtual knot from its gauss code.

Consists of two sub-scripts:
- `src/dinv` calculates a d-invariant from a tait graph.
- `src/tait` calculates a tait graph from a gauss code.

# Usage

- Compute d-invariant of some knots:

		python main.py <file>

## Input
Text file with each line of format `<name> <gauss code>`.

where:

- `<name>` is any string containing no whitespace.
- `<gauss code>` is a dash-delimited sequence describing the number of crossing traversed, whether traversed over/under and sign of the crossing when traversing the knot.

Example is included in `/example_data/knots.txt`.

## Output

The `/out` directory is cleared (at surface level) and for each knot three files are created:

- `<name>.tait`

Edge list of each tait graph spearated by empty line.


- `<name>.dinv`

Sympy unicode output of map of d-invariant of lattice integer flows for each tait graph, separated by empty line.


- `<name>.dinv.pickle`

Pickled tuple of the two maps above.

# Dependencies

The following python packages:
- `networkx v2.6.3`
- `progress v1.6`
- `sympy v1.11`

#!/bin/bash
cd "$(dirname $0)"

# induced_matrix_and_tree
rm -f data/pruned-A-Dultrametric.tre data/pruned-A-Daminoacid.fas
dendrobites/induced_matrix_and_tree.py --char=data/A-Daminoacid.fas --tree=data/A-Dultrametric.tre  --data-type=protein A B C || exit
diff data/pruned-A-Dultrametric.tre test/expected/pruned-A-Dultrametric.tre || exit
diff data/pruned-A-Daminoacid.fas test/expected/pruned-A-Daminoacid.fas || exit
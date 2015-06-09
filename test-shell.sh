#!/bin/bash
cd "$(dirname $0)"

# induced_matrix_and_tree
rm -f data/pruned-A-Dultrametric.tre data/pruned-A-Daminoacid.fas
dendrobites/induced_matrix_and_tree.py --char=data/A-Daminoacid.fas --tree=data/A-Dultrametric.tre  --data-type=protein A B C || exit
diff data/pruned-A-Dultrametric.tre test/expected/pruned-A-Dultrametric.tre || exit
diff data/pruned-A-Daminoacid.fas test/expected/pruned-A-Daminoacid.fas || exit

# induced_matrix_and_tree
rm -f test/output/tip-match-correct test/output/tip-match-error
python dendrobites/tip_label_match.py --char data/A-Dnucleotide.fas --tree data/A-Dultrametric.tre > test/output/tip-match-correct || exit
python dendrobites/tip_label_match.py --char data/A-Dnucleotide_label_error.fas --tree data/A-Dultrametric.tre 2> test/output/tip-match-error || exit
diff test/output/tip-match-correct  test/expected/tip-match-correct || exit
diff test/output/tip-match-error test/expected/tip-match-error || exit
#!/usr/bin/env python
'''Read in a matrix and a proportion of invariant sites.
For each pair of taxa calculate a pairwise alignment length
    (this is simply the total alignment length minus the
     number of columns for which both members of the pair
     have gaps).
Use the average pairwise alignment length as an estimator
    of the equilibrium length of sequences evolving under the
    paired-invariants model of McTavish, Steel, Holder
    http://de.arxiv.org/abs/1504.07124

Given the proportion of invariant sites (command line argument)
    and this estimate of the of equilibrium length, estimate
    the proportion of gapless, constant columns that are in
    the invariant class. Call this `p_r`

Finally, walk through the matrix and identify all of the constant
    gapless patterns. For each one, retain it with a probability of
    1 - p_r

Write the resulting matrix out. If the data were simulated under
    the paired-invariants model, and the sequences are long enough,
    then the equilibrium length should be estimated accurately.
    If the command-line argument for the generating p_invariant is
    correct, then p_r should be accurately estimated. In this case (of
    no model misspecification), MTH conjectures that the resulting
    matrix should be yield a consistent matrix of the tree under a
    distance correction that uses no rate heterogeneity.
'''
from dendropy.datamodel.charmatrixmodel import data_type_matrix_map, \
                                               DnaCharacterMatrix


def induced_matrix_and_tree(char_mat_filepath,
                            tree_filepath,
                            taxa_labels,
                            char_type=DnaCharacterMatrix,
                            char_schema='fasta',
                            tree_schema='newick'):
    '''Reads an (optional) CharacterMatrix from `char_mat_filepath` and
    a (required) tree from `tree_filepath`. Prunes both down to just
    the taxa whose labels match `taxa_labels` and then returns (char_mat, tree).
    '''
    # read the char matrix and tree....
    char_mat, tree = read_matrix_and_tree(char_mat_filepath,
                                          tree_filepath,
                                          char_type=char_type,
                                          char_schema=char_schema,
                                          tree_schema=tree_schema)
    taxa_labels = frozenset(taxa_labels)
    if not tree.taxon_namespace.has_taxa_labels(taxa_labels):
        for t in taxa_labels:
            if not tree.taxon_namespace.has_taxon_label(t):
                raise ValueError('Taxon "{}" not found in the taxon namespace of this data.\n'.format(t))
    to_cull = []
    for t in tree.taxon_namespace:
        if t.label not in taxa_labels:
            to_cull.append(t)
    if to_cull:
        tree.prune_taxa(to_cull)
        if char_mat:
            char_mat.remove_sequences(to_cull)
    return char_mat, tree

def get_path_with_prefix(template, filename_prefix):
    directory, fn = os.path.split(os.path.abspath(template))
    ocfn = filename_prefix + fn
    return os.path.join(directory, ocfn)

def _main(char_mat_filepath,
          data_type_name,
          p_inv,
          char_schema='nexus'):
    # Validate the data_type argument and use it to find the CharacterMatrix type
    dt = data_type_name.lower()
    mat_type = data_type_matrix_map.get(dt)
    if mat_type is None:
        emf = 'The data type "{u}" is not recognized.\nExpecting one of "{t}".\n'
        k = data_type_matrix_map.keys()
        k.sort()
        raise ValueError(emf.format(u=data_type_name, t='", "'.join(k)))
    # read the char matrix 
    char_mat = mat_type.get(path=char_mat_filepath, schema=char_schema)


if __name__ == '__main__':
    import argparse
    import sys
    import os
    script_name = os.path.split(sys.argv[0])[1]
    description = '''Subsample constant, gapless columns as if they were generated under the paired-invariants model.'''
    parser = argparse.ArgumentParser(prog=script_name, description=description)
    parser.add_argument('--data-type', default='dna', type=str, required=False, help='a data_type. Default is "dna"')
    parser.add_argument('datafile', default=None, nargs=1, help='filepath of the character data')
    parser.add_argument('--p-inv', required=True, type=float, help='A proportion of invariant sites for the paired-invariants model')
    args = parser.parse_args(sys.argv[1:])
    try:
        assert args.p_inv > 0.0
        assert args.p_inv < 1.0
        assert len(args.datafile) == 1
        _main(args.datafile[0], args.data_type, p_inv=args.p_inv)
    except Exception as x:
        raise
        sys.exit('{}: {}\n'.format(script_name, str(x)))

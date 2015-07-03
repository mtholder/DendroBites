#!/usr/bin/env python
'''Read in a matrix and a proportion of invariant sites.

Uses the average pairwise sequence length as an estimator
    of the equilibrium length of sequences evolving under the
    paired-invariants model of McTavish, Steel, Holder
    http://de.arxiv.org/abs/1504.07124

Given the proportion of invariant sites (command line argument)
    and this estimate of the of equilibrium length, estimate
    the proportion of gapless, constant columns that are in
    the invariant class. Call this `p_r`

Retain approximately (1 - p_r) of constant, gapless columns. This
    is approximate because, if the number of input constant, gapless
    columns is C then the proportions attainable exactly by deletion
    are 1/C, 2/C, ... (C-1)/C, 1.0

The script attempts to cull the columns in proportion to the frequency
    of the states among the constant gapless columns so that it
    will distort the state frequencies of these columns as little as
    possible.

Write the resulting matrix out in the same schema as the input data.

If :
  1. the data were simulated under the paired-invariants model, 
  2. the sequences are long enough, and
  3. the --p-inv argument accurately provides the proportion of
    invariant sites
then
  1. the equilibrium length should be estimated accurately,
  2. the matrix written out should be a reasonable proxy for the
    free-to-vary subset of the matrix
  3. MTH conjectures that preprocessing of matrices and analyzing the
    resulting matrices should provide a consistent estimate of
    phylogeny and branch lengths (provided that the substitution
    model and original tree correspond to a model that yields 
    consistent estimates).
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

def iter_columns(char_mat, taxa_order=None):
    '''Iterates through the columns of a `char_mat`. Returning 
    the cells in a list with the order determined by
    `taxa_order` iterable (or the order of taxa in the char_mat if `None`)
    '''
    if taxa_order is None:
        taxa_order = [i for i in char_mat]
    full_rows = [char_mat[i] for i in taxa_order]
    nc = len(full_rows[0])
    for row in full_rows:
        if len(row) != nc:
            raise ValueError('iter_columns requires aligned matrices.')
    for i in xrange(nc):
        column = []
        for row in full_rows:
            column.append(row[i])
        yield column

def is_gapless_constant(column):
    '''Takes and iterable of States. Returns `True` if there are
    no gaps, and the same symbol is used by every element in `column`.
    '''
    first_sym = None
    for cell in column:
        if cell.is_gap_state:
            return False
        if first_sym is None:
            first_sym = cell.symbol
        elif cell.symbol != first_sym:
            return False
    return True

def count_gap_cells(column):
    '''Returns the number of cells in `column` for which is_gap_state is `True`'''
    x = 0
    for c in column:
        if c.is_gap_state:
            x += 1
    return x

def characterize_mat_wrt_const_gapless(char_mat):
    '''Walks through `char_mat`
    returns:
       1. the total # of columns,
       2. the number of cells that are gaps
       3. a map of a state symbol to the set of column indices for
        the columns that are constant (and gapless) for that symbol.
    '''
    n_col = 0
    const_col_type2ind_set = {}
    num_gap_cells = 0
    for col_ind, column in enumerate(iter_columns(char_mat)):
        if is_gapless_constant(column):
            symbol = column[0].symbol
            ind_set = const_col_type2ind_set.get(symbol)
            if ind_set is None:
                ind_set = set()
                const_col_type2ind_set[symbol] = ind_set
            ind_set.add(col_ind)
        else:
            num_gap_cells += count_gap_cells(column)
        n_col = 1 + col_ind
    return (n_col, num_gap_cells, const_col_type2ind_set)

def calc_num_to_cull_by_state(num_inv_columns, symbol2ind_set):
    '''Takes `symbol2ind_set` which maps state symbols to iterable
    collections of indices that are a partition of the full set of indices.

    Returns a mapping of these symbols to a pair of
        1. the # of items to remove from this state
        2. the total number of items for this state on input.

    The numbers to remove for each set of indices are chosen such that:
            1. the total number of times removed is as close to `num_inv_column`
                as feasible, and 
            2. the # of indices removed for a state is proportional to relative
                size of the input index set for that state (the proportions of the
                the states post-pruning will be close to the original proportion of
                symbols).
    '''
    num_const_gapless = sum([len(v) for v in symbol2ind_set.values()])
    invariant_frac = num_inv_columns/float(num_const_gapless)
    ideal_num_to_cull = int(round(num_inv_columns))
    num_to_cull_by_state = {}
    num_left_to_cull = ideal_num_to_cull
    for state, col_ind_for_this_state in symbol2ind_set.items():
        num_to_cull_for_this_state = int(round(invariant_frac*len(col_ind_for_this_state)))
        num_to_cull_for_this_state = min(len(col_ind_for_this_state), num_to_cull_for_this_state)
        num_to_cull_for_this_state = min(num_to_cull_for_this_state, num_left_to_cull)
        num_left_to_cull -= num_to_cull_for_this_state
        num_to_cull_by_state[state] = (num_to_cull_for_this_state, len(col_ind_for_this_state))
    # Deal with rounding error
    if num_left_to_cull > 0:
        sym_list = num_to_cull_by_state.keys()
        sym_list.sort()
        for state in sym_list:
            tc, tot = num_to_cull_by_state[state]
            ntc = min(tc + num_left_to_cull, tot)
            if ntc != tc:
                num_added = ntc - tc
                num_to_cull_by_state[state] = (ntc, tc)
                num_left_to_cull -= num_added
                if num_left_to_cull == 0:
                    break
    return num_to_cull_by_state

def create_inds_to_cull_from_numbers_to_cull(symbol2ind_set, num_to_cull_by_state):
    total_to_cull = set()
    for state, col_ind_for_this_state in symbol2ind_set.items():
        num_to_cull_for_this_state = num_to_cull_by_state[state][0]
        to_cull = set()
        for ind in col_ind_for_this_state:
            if len(to_cull) < num_to_cull_for_this_state:
                to_cull.add(ind)
            else:
                break
        total_to_cull.update(to_cull)
    return total_to_cull

def calc_inds_to_cull(num_inv_columns, const_col_type2ind_set):
    num_to_cull_by_state = calc_num_to_cull_by_state(num_inv_columns=num_inv_columns,
                                                     symbol2ind_set=const_col_type2ind_set)
    return create_inds_to_cull_from_numbers_to_cull(symbol2ind_set=const_col_type2ind_set,
                                                    num_to_cull_by_state=num_to_cull_by_state)

def new_mat_by_del_paired_invariants(char_mat, p_inv):
    '''Takes a char_mat that is assumed to be a product of evolution by the paired-invariants
    model with a proportion of invariant sites equal to p_inv.
    Returns a proxy for the a matrix representing the results of the free-to-vary evolution
    by removing constant, gapless columns from char_mat.
    '''
    r = characterize_mat_wrt_const_gapless(char_mat)
    num_cols, num_gap_cells, const_col_type2ind_set = r
    num_taxa = len(char_mat)
    est_equil_len = (num_cols*len(char_mat) - num_gap_cells)/float(num_taxa)
    est_num_inv_columns = p_inv*est_equil_len
    to_cull = calc_inds_to_cull(num_inv_columns=est_num_inv_columns,
                                const_col_type2ind_set=const_col_type2ind_set)
    to_retain = set(xrange(num_cols))
    to_retain.difference_update(to_cull)
    label = 'to_retain'
    char_mat.new_character_subset(label=label, character_indices=to_retain)
    retained = char_mat.export_character_subset(character_subset=label)
    # need to remove the char subset because they are not being re-indexed.
    del retained.character_subsets[label]
    return retained


def _main(char_mat_filepath,
          data_type_name,
          p_inv,
          schema='nexus'):
    # Validate the data_type argument and use it to find the CharacterMatrix type
    dt = data_type_name.lower()
    mat_type = data_type_matrix_map.get(dt)
    if mat_type is None:
        emf = 'The data type "{u}" is not recognized.\nExpecting one of "{t}".\n'
        k = data_type_matrix_map.keys()
        k.sort()
        raise ValueError(emf.format(u=data_type_name, t='", "'.join(k)))
    # read the char matrix 
    char_mat = mat_type.get(path=char_mat_filepath, schema=schema)
    retained = new_mat_by_del_paired_invariants(char_mat, p_inv)
    retained.write_to_stream(sys.stdout, schema=schema)

if __name__ == '__main__':
    import argparse
    import sys
    import os
    script_name = os.path.split(sys.argv[0])[1]
    description = '''Subsample constant, gapless columns as if they were generated under the paired-invariants model.'''
    parser = argparse.ArgumentParser(prog=script_name, description=description)
    parser.add_argument('--data-type', default='dna', type=str, required=False, help='a data_type. Default is "dna"')
    parser.add_argument('--schema', default='nexus', type=str, required=False, help='A file format name. Default is "nexus"')
    parser.add_argument('datafile', default=None, nargs=1, help='filepath of the character data')
    parser.add_argument('--p-inv', required=True, type=float, help='A proportion of invariant sites for the paired-invariants model')
    args = parser.parse_args(sys.argv[1:])
    try:
        assert args.p_inv > 0.0
        assert args.p_inv < 1.0
        assert len(args.datafile) == 1
        _main(args.datafile[0], args.data_type, schema=args.schema, p_inv=args.p_inv)
    except Exception as x:
        raise
        sys.exit('{}: {}\n'.format(script_name, str(x)))

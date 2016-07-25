#!/usr/bin/env python
'''Read in a matrix and a list of taxon labels.

Prints out all of the columns of the character matrix for which
the indicated taxa have a set of states that does not overlap 
with the symbols displayed by the "outgroup" (the set of 
taxa that are not listed).

This identifies sites that are potentially clean synapomorphies
for the listed taxa.

Cells with missing data in the matrix are ignored. 
'''
from dendropy.datamodel.charmatrixmodel import data_type_matrix_map, \
                                               DnaCharacterMatrix

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

def find_potential_synapo_columns(char_mat, ingroup_taxa):
    r = []
    all_taxa = char_mat.taxon_namespace
    for col_ind, column in enumerate(iter_columns(char_mat)):
        in_c = set()
        out_c = set()
        are_disjunct = True
        for tax_ind, cell in enumerate(column):
            if cell.is_gap_state or (not cell.is_single_state):
                continue
            symbol = cell.symbol
            taxon = all_taxa[tax_ind]
            if taxon in ingroup_taxa:
                if symbol in out_c:
                    are_disjunct = False
                    break
                else:
                    in_c.add(symbol)
            else:
                if symbol in in_c:
                    are_disjunct = False
                    break
                else:
                    out_c.add(symbol)
        if are_disjunct and bool(in_c) and bool(out_c):
            r.append([col_ind, in_c, out_c])
    return r

def _main(char_mat_filepath,
          data_type_name,
          taxa_identifiers,
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
    tns = char_mat.taxon_namespace
    nt = len(tns)
    rts = set(taxa_identifiers)
    if len(rts) != len(taxa_identifiers):
        raise ValueError("Some taxa labels were repeated")
    ingroup_taxa = frozenset(tns.get_taxa(taxa_identifiers, first_match_only=True))
    if len(ingroup_taxa) != len(taxa_identifiers):
        fts = {i.label for i in ingroup_taxa}
        msg = '", "'.join([i for i in rts - fts])
        raise ValueError('Could not find the taxa labels: "{}"'.format())
    if len(ingroup_taxa) == nt:
        raise ValueError('Listing all tips is nonsensical')
    psc = find_potential_synapo_columns(char_mat, ingroup_taxa)
    for el in psc:
        col_in, in_states, out_states = el
        sys.stdout.write('Column {}: in states = {{{}}}. out states = {{{}}}.\n'.format(col_in,
                                                                                ', '.join([i for i in in_states]),
                                                                                ', '.join([i for i in out_states])))
if __name__ == '__main__':
    import argparse
    import sys
    import os
    script_name = os.path.split(sys.argv[0])[1]
    description = '''Find sites that are putative synapomorphies for the taxa indicated.'''
    parser = argparse.ArgumentParser(prog=script_name, description=description)
    parser.add_argument('--data-type', default='dna', type=str, required=False, help='a data_type. Default is "dna"')
    parser.add_argument('--char-mat', type=str, required=True, help='A filepath for the input file')
    parser.add_argument('--schema', default='nexus', type=str, required=False, help='A file format name. Default is "nexus"')
    parser.add_argument('taxa', default=None, nargs='+', help='list of taxon names for the group whose synapomorphies that you want to find')
    args = parser.parse_args(sys.argv[1:])
    try:
        assert len(args.taxa) > 0
        _main(args.char_mat, args.data_type, taxa_identifiers=args.taxa, schema=args.schema)
    except Exception as x:
        raise
        sys.exit('{}: {}\n'.format(script_name, str(x)))

#!/usr/bin/env python
'''Using dendropy to prune a tree to an induced tree and
prune the same set of removed taxa from a data matrix.
'''
from dendropy import Tree
from dendropy.datamodel.charmatrixmodel import data_type_matrix_map, \
                                               DnaCharacterMatrix
def read_matrix_and_tree(char_file_path,
                         tree_file_path,
                         char_type=DnaCharacterMatrix,
                         char_schema='fasta',
                         tree_schema='newick'):
    if char_file_path:
        d = char_type.get(path=char_file_path, schema=char_schema)
        tn = d.taxon_namespace
        tn.is_mutable = False
    else:
        d, tn = None, None
    tree = Tree.get(path=tree_file_path,
                    schema=tree_schema,
                    preserve_underscores=True,
                    taxon_namespace=tn)
    treed_taxa = [i.taxon for i in tree.leaf_nodes()]
    if len(treed_taxa) != len(d.taxon_namespace):
        missing = [i.label for i in d.taxon_namespace if i not in treed_taxa]
        emf = 'Some of the taxa are not in the tree. Missing "{}"\n'
        em = emf.format('", "'.join(missing))
        raise ValueError(em)
    return d, tree

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
          tree_filepath,
          taxa_labels,
          data_type_name,
          char_schema='fasta',
          tree_schema='newick'):
    # Validate the data_type argument and use it to find the CharacterMatrix type
    dt = data_type_name.lower()
    mat_type = data_type_matrix_map.get(dt)
    if mat_type is None:
        emf = 'The data type "{u}" is not recognized.\nExpecting one of "{t}".\n'
        k = data_type_matrix_map.keys()
        k.sort()
        raise ValueError(emf.format(u=data_type_name, t='", "'.join(k)))
    # Validate the output filenames and make sure we won't overwrite content.
    out_paths = []
    if char_mat_filepath:
        out_char = get_path_with_prefix(char_mat_filepath, 'pruned-')
        out_paths.append(out_char)
    out_tree = get_path_with_prefix(tree_filepath, 'pruned-')
    out_paths.append(out_tree)
    for ofp in out_paths:
        if os.path.exists(ofp):
            raise RuntimeError('"{}" already exists! Move it before running this script.\n'.format(ofp))
    # read the char matrix and tree....
    char_mat, tree = induced_matrix_and_tree(char_mat_filepath,
                                             tree_filepath,
                                             taxa_labels,
                                             char_type=mat_type,
                                             char_schema=char_schema,
                                             tree_schema=tree_schema)
    tree.write_to_path(out_tree, schema=tree_schema)
    if char_mat:
        char_mat.write_to_path(out_char, schema=char_schema)

if __name__ == '__main__':
    import argparse
    import sys
    import os
    script_name = os.path.split(sys.argv[0])[1]
    description = '''Takes a data file and a tree and a series of taxon labels.
Writes pruned versions of the matrix and tree out to files with a "pruned-" prefix
and the same ending as the input files.'''
    parser = argparse.ArgumentParser(prog=script_name, description=description)
    parser.add_argument('--data-type', default='dna', type=str, required=False, help='a data_type. Default is "dna"')
    parser.add_argument('--char', default=None, type=str, required=False, help='filepath of the character data')
    parser.add_argument('--tree', default=None, type=str, required=Tree, help='filepath of the tree')
    parser.add_argument('taxa', nargs='+')
    args = parser.parse_args(sys.argv[1:])
    try:
        _main(args.char, args.tree, args.taxa, args.data_type)
    except Exception as x:
        sys.exit('{}: {}\n'.format(script_name, str(x)))

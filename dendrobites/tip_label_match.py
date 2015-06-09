#!/usr/bin/env python
'''Check that taxa in alignment match those in tree,
if not report differences'''
from dendropy import DnaCharacterMatrix, Tree
from dendropy.datamodel.charmatrixmodel import data_type_matrix_map, \
                                               DnaCharacterMatrix

def mutable_read_matrix_and_tree(char_file_path,
                                tree_file_path,
                                char_type=DnaCharacterMatrix,
                                char_schema='fasta',
                                tree_schema='newick'):
    '''Reads in tree and character matrix,
    mutable namespace means names may not match'''
    if char_file_path:
        char_mat = char_type.get(path=char_file_path,
                               schema=char_schema)
        # make the taxon_namespace mutable,
        # so that tree can be read even if different
        char_mat.taxon_namespace.is_mutable = True
        tree = Tree.get(path=tree_file_path,
                        schema=tree_schema,
                        preserve_underscores=True,
                        taxon_namespace=char_mat.taxon_namespace)
    else:
        char_mat, tree = None, None
    return char_mat, tree

def tip_label_match(char_mat_filepath,
                            tree_filepath,
                            char_type=DnaCharacterMatrix,
                            char_schema='fasta',
                            tree_schema='newick'):
    '''Reads a (required) CharacterMatrix from `char_mat_filepath`
    and a (required) tree from `tree_filepath` and
    checks if tip labels match'''
    char_mat, tree = mutable_read_matrix_and_tree(char_mat_filepath,
                                          tree_filepath,
                                          char_type=DnaCharacterMatrix,
                                          char_schema=char_schema,
                                          tree_schema=tree_schema)
    treed_taxa = [i.taxon for i in tree.leaf_nodes()]
    if set(treed_taxa) != char_mat.poll_taxa():
        tree_missing = [i.label for i in char_mat.taxon_namespace if i not in treed_taxa]
        emf = 'Some of the taxa in the matrix are not in the tree.\
                Tree is missing "{}"\n'
        em = emf.format('", "'.join(tree_missing))
#        raise ValueError(em)
        sys.stderr.write(em)
        mat_missing = [i.label for i in char_mat.taxon_namespace if i not in char_mat.poll_taxa()]
        emf = 'Some of the taxa in the tree are not in the data matrix.\
                Matrix is missing "{}"\n'
        em = emf.format('", "'.join(mat_missing))
 #       raise ValueError(em)
        sys.stderr.write(em)
        return 0
    else:
        return 1



def _main(char_mat_filepath,
          tree_filepath,
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
    # read the char matrix and tree....
    match_check = tip_label_match(char_mat_filepath,
                    tree_filepath,
                    char_type=data_type_name,
                    char_schema=char_schema,
                    tree_schema=tree_schema)
    if match_check:
        sys.stdout.write("Tips match\n")



if __name__ == '__main__':
    import argparse
    import sys
    import os
    script_name = os.path.split(sys.argv[0])[1]
    description = '''Takes a data file and a tree.
    If taxon labels aren't matched, returns labels that are found in only one or the other.'''
    parser = argparse.ArgumentParser(prog=script_name, description=description)
    parser.add_argument('--char', default=None, type=str, required=True, help='filepath of the character data')
    parser.add_argument('--tree', default=None, type=str, required=True, help='filepath of the tree')
    parser.add_argument('--char-schema', default="fasta", type=str, required=False, help='schema for the character data')
    parser.add_argument('--tree-schema', default="newick", type=str, required=False, help='schema for the tree')
    parser.add_argument('--data-type', default='dna', type=str, required=False, help='a data_type. Default is "dna"') #This doesn't actaully do anything...
    args = parser.parse_args(sys.argv[1:])
    try:
        _main(args.char, args.tree, args.data_type, args.char_schema, args.tree_schema)
    except Exception as x:
        sys.exit('{}: {}\n'.format(script_name, str(x)))

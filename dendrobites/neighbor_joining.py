#!/usr/bin/env python
import dendropy
from dendropy.calculate.phylogeneticdistance import PhylogeneticDistanceMatrix as DendropyDistMat
def parse_distances(fn):
    mat = {}
    with open(fn, 'rU') as inp:
        for line in inp:
            bogus, first, second, n, dist = line.strip().split()
            f, s = int(first), int(second)
            d = float(dist)
            mat.setdefault(f, {})[s] = d
            mat.setdefault(s, {})[f] = d
    return mat
def _main(jkk_ssv_filepath):
    # parse to a dict of dicts with integer "names" as the keys
    dist_mat = parse_distances(jkk_ssv_filepath)
    # Convert it to a special PhylogeneticDistanceMatrix from dendropy
    dendropy_dist = convert_to_dendropy_dist(dist_mat)
    nj = dendropy_dist.nj_tree(is_weighted_edge_distances=True)
    nj.print_plot(plot_metric='length')

def convert_to_dendropy_dist(dist_mat):
    '''Takes a distance matrix as a dict of dicts.
    Creates a taxon namespace for the keys and then creates a PhylogeneticDistanceMatrix
    from the distances
    '''
    taxon_namespace = dendropy.TaxonNamespace(label="taxa")
    names = list(dist_mat.keys())
    names.sort()
    name2taxon = {}
    for taxname in names:
        taxon = dendropy.Taxon(str(taxname))
        name2taxon[taxname] = taxon
        taxon_namespace.add_taxon(taxon)
    by_taxa = {}
    for tax_name in names:
         taxon = name2taxon[tax_name]
         btr = {taxon: 0.0}
         by_taxa[taxon] = btr
         row = dist_mat[tax_name]
         for other_name, d in row.items():
            other_taxon = name2taxon[other_name]
            btr[other_taxon] = d   
    
    dendropy_dist = DendropyDistMat()
    dendropy_dist.compile_from_dict(by_taxa, taxon_namespace)
    return dendropy_dist

if __name__ == '__main__':
    import argparse
    import sys
    import os
    script_name = os.path.split(sys.argv[0])[1]
    description = '''Takes a filepath to a quirky space-separated representation of the distance matrix. Prints an NJ tree.'''
    parser = argparse.ArgumentParser(prog=script_name, description=description)
    parser.add_argument('distances')
    args = parser.parse_args(sys.argv[1:])
    try:
        _main(args.distances)
    except Exception as x:
        sys.exit('{}: {}\n'.format(script_name, str(x)))

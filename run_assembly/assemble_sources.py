import sys
import random
from util import *
import indra.tools.assemble_corpus as ac

if __name__ == '__main__':
    #sources = ['bel', 'biopax', 'reach', 'sparser']
    cmd = sys.argv[1]
    # Break each assembly cmd down into intermediate steps
    # Step 1: Load the statement files and combine into a single pickle
    if cmd == 'load_stmts':
        sources = [
            'output/bioexp_bel.pkl',
            'output/bioexp_biopax_fixed.pkl',
            'data/bioexp_sparser.pkl']
        stmts = []
        for source in sources:
            stmts += ac.load_statements(source)
        pkldump(stmts, 'all_raw') # FIXME
    # FIXME: Remove
    elif cmd == 'sample':
        stmts = pklload('all_raw')
        random.shuffle(stmts)
        stmts = stmts[0:10000]
        pkldump(stmts, 'all_raw_sample')
    elif cmd == 'filter_no_hypothesis':
        stmts = pklload('all_raw_sample') # FIXME
        stmts = ac.filter_no_hypothesis(stmts, save=prefixed_pkl(cmd))
    elif cmd == 'map_grounding':
        stmts = pklload('filter_no_hypothesis')
        stmts = ac.map_grounding(stmts, save=prefixed_pkl(cmd))
    elif cmd == 'filter_grounded_only':
        stmts = pklload('map_grounding')
        stmts = ac.filter_grounded_only(stmts, save=prefixed_pkl(cmd))
    elif cmd == 'filter_genes_only':
        stmts = pklload('filter_grounded_only')
        stmts = ac.filter_genes_only(stmts, specific_only=False,
                                     save=prefixed_pkl(cmd))
    elif cmd == 'filter_human_only':
        stmts = pklload('filter_genes_only')
        stmts = ac.filter_human_only(stmts, save=prefixed_pkl(cmd))
    elif cmd == 'expand_families':
        stmts = pklload('filter_human_only')
        stmts = ac.expand_families(stmts, save=prefixed_pkl(cmd))
    elif cmd == 'filter_gene_list':
        stmts = pklload('expand_families')
        # Load the genes from the outer build directory
        with open('../build/prior_genes.txt', 'rt') as f:
            data_genes = [line.strip() for line in f.readlines()]
        stmts = ac.filter_gene_list(stmts, data_genes, 'one',
                                    save=prefixed_pkl(cmd))
    elif cmd == 'map_sequence':
        stmts = pklload('filter_gene_list')
        stmts = ac.map_sequence(stmts, save=prefixed_pkl(cmd))
    elif cmd == 'run_preassembly':
        stmts = pklload('map_sequence')
        # TODO: Set poolsize for full assembly run
        stmts = ac.run_preassembly(stmts, return_toplevel=False,
                                   save=prefixed_pkl(cmd),
                                   poolsize=16)
    """
    stmts = ac.filter_no_hypothesis(stmts)
    stmts = ac.map_grounding(stmts, save=prefixed_pkl('grounded'))
    stmts = ac.filter_grounded_only(stmts)
    stmts = ac.filter_genes_only(stmts, specific_only=False)    
    stmts = ac.filter_human_only(stmts)
    stmts = ac.expand_families(stmts)
    stmts = ac.filter_gene_list(stmts, data_genes, 'one')
    stmts = ac.map_sequence(stmts, save=pjoin(outf, 'smapped.pkl'))
    stmts = ac.run_preassembly(stmts, return_toplevel=False,
                               save=pjoin(outf, 'preassembled.pkl'),
                               poolsize=12)
    """

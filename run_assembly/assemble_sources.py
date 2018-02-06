import sys
import random
from util import *
import indra.tools.assemble_corpus as ac
from indra.sources import signor
from indra.mechlinker import MechLinker

if __name__ == '__main__':
    #sources = ['bel', 'biopax', 'reach', 'sparser']
    cmd = sys.argv[1]
    # Break each assembly cmd down into intermediate steps
    if cmd == 'signor':
        sp = signor.SignorProcessor()
        pkldump(sp.statements, 'signor')
    elif cmd == 'load_stmts':
        sources = [
            'output/bioexp_signor.pkl',
            'output/bioexp_bel.pkl',
            'output/bioexp_biopax_fixed.pkl',
            'data/bioexp_reach.pkl',
            'data/bioexp_sparser.pkl']
        stmts = []
        for source in sources:
            stmts += ac.load_statements(source)
        pkldump(stmts, 'all_raw')
    elif cmd == 'filter_no_hypothesis':
        stmts = pklload('all_raw')
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
    elif cmd == 'preassembled':
        stmts = pklload('map_sequence')
        stmts = ac.run_preassembly(stmts, return_toplevel=False,
                                   save=prefixed_pkl(cmd),
                                   poolsize=16)
    elif cmd == 'filter_belief':
        stmts = pklload('preassembled')
        stmts = ac.filter_belief(stmts, 0.95, save=prefixed_pkl(cmd))
    elif cmd == 'filter_top_level':
        stmts = pklload('filter_belief')
        stmts = ac.filter_top_level(stmts, save=prefixed_pkl(cmd))
    elif cmd == 'filter_enzyme_kinase':
        stmts = pklload('filter_top_level')
        stmts = ac.filter_enzyme_kinase(stmts, save=prefixed_pkl(cmd))
    elif cmd == 'filter_mod_nokinase':
        stmts = pklload('filter_enzyme_kinase')
        stmts = ac.filter_mod_nokinase(stmts, save=prefixed_pkl(cmd))
    elif cmd == 'reduce_activities':
        stmts = pklload('filter_mod_nokinase')
        ml = MechLinker(stmts)
        ml.gather_explicit_activities()
        ml.reduce_activities()
        pkldump(stmts, cmd)
    elif cmd == 'reduce_mods':
        stmts = pklload('reduce_activities')
        ml = MechLinker(stmts)
        ml.gather_modifications()
        ml.reduce_modifications()
        pkldump(stmts, cmd)


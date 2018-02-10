import sys
import random
from util import *
import indra.tools.assemble_corpus as ac
from indra.sources import signor
from indra.mechlinker import MechLinker

if __name__ == '__main__':
    #sources = ['bel', 'biopax', 'reach', 'sparser']
    cmd = sys.argv[1]
    assembly_cmds = (
            'filter_no_hypothesis', 'map_grounding', 'filter_grounded_only',
            'filter_genes_only', 'filter_human_only', 'expand_families',
            'filter_gene_list', 'map_sequence', 'preassembled',
            'filter_belief', 'filter_top_level', 'filter_enzyme_kinase',
            'filter_mod_nokinase', 'reduce_activities', 'reduce_mods')
    if cmd in assembly_cmds:
        input_file = sys.argv[2]
        output_file = sys.argv[3]
    # Break each assembly cmd down into intermediate steps
    if cmd == 'signor':
        sp = signor.SignorProcessor()
        pkldump(sp.statements, 'signor')
    elif cmd == 'load_stmts':
        sources = [
            'output/bioexp_signor.pkl',
            'output/bioexp_bel.pkl',
            'data/bioexp_biopax_fixed.pkl',
            'data/bioexp_reach.pkl',
            'data/bioexp_sparser.pkl']
        stmts = []
        for source in sources:
            stmts += ac.load_statements(source)
        pkldump(stmts, 'all_raw')
    elif cmd == 'filter_no_hypothesis':
        stmts = pklload(input_file)
        stmts = ac.filter_no_hypothesis(stmts, save=prefixed_pkl(output_file))
    elif cmd == 'map_grounding':
        stmts = pklload(input_file)
        stmts = ac.map_grounding(stmts, save=prefixed_pkl(output_file))
    elif cmd == 'filter_grounded_only':
        stmts = pklload(input_file)
        stmts = ac.filter_grounded_only(stmts, save=prefixed_pkl(output_file))
    elif cmd == 'filter_genes_only':
        stmts = pklload(input_file)
        stmts = ac.filter_genes_only(stmts, specific_only=False,
                                     save=prefixed_pkl(output_file))
    elif cmd == 'filter_human_only':
        stmts = pklload(input_file)
        stmts = ac.filter_human_only(stmts, save=prefixed_pkl(output_file))
    elif cmd == 'expand_families':
        stmts = pklload(input_file)
        stmts = ac.expand_families(stmts, save=prefixed_pkl(output_file))
    elif cmd == 'filter_gene_list':
        stmts = pklload(input_file)
        # Load the genes from the outer build directory
        with open('../build/prior_genes.txt', 'rt') as f:
            data_genes = [line.strip() for line in f.readlines()]
        stmts = ac.filter_gene_list(stmts, data_genes, 'one',
                                    save=prefixed_pkl(output_file))
    elif cmd == 'map_sequence':
        stmts = pklload(input_file)
        stmts = ac.map_sequence(stmts, save=prefixed_pkl(output_file))
    elif cmd == 'preassembled':
        stmts = pklload(input_file)
        stmts = ac.run_preassembly(stmts, return_toplevel=False,
                                   save=prefixed_pkl(output_file),
                                   poolsize=16)
    elif cmd == 'filter_belief':
        stmts = pklload(input_file)
        stmts = ac.filter_belief(stmts, 0.95, save=prefixed_pkl(output_file))
    elif cmd == 'filter_top_level':
        stmts = pklload(input_file)
        stmts = ac.filter_top_level(stmts, save=prefixed_pkl(output_file))
    elif cmd == 'filter_enzyme_kinase':
        stmts = pklload(input_file)
        stmts = ac.filter_enzyme_kinase(stmts, save=prefixed_pkl(output_file))
    elif cmd == 'filter_mod_nokinase':
        stmts = pklload(input_file)
        stmts = ac.filter_mod_nokinase(stmts, save=prefixed_pkl(output_file))
    elif cmd == 'reduce_activities':
        stmts = pklload(input_file)
        ml = MechLinker(stmts)
        ml.gather_explicit_activities()
        ml.reduce_activities()
        pkldump(ml.statements, output_file)
    elif cmd == 'reduce_mods':
        stmts = pklload(input_file)
        ml = MechLinker(stmts)
        ml.gather_modifications()
        ml.reduce_modifications()
        pkldump(ml.statements, output_file)


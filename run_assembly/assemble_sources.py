import sys
from util import *
import indra.tools.assemble_corpus as ac
from indra.sources import signor
from indra.mechlinker import MechLinker
from scorer import CuratedScorer

if __name__ == '__main__':
    cmd = sys.argv[1]
    assembly_cmds = (
            'filter_no_hypothesis', 'map_grounding', 'filter_grounded_only',
            'filter_genes_only', 'filter_human_only', 'expand_families',
            'map_sequence', 'preassembled')
    if cmd in assembly_cmds:
        input_file = sys.argv[2]
        output_file = sys.argv[3]
    # Break each assembly cmd down into intermediate steps
    if cmd == 'signor':
        sp = signor.process_from_web()
        pkldump(sp.statements, 'signor')
    elif cmd == 'load_stmts':
        output = sys.argv[2]
        sources = sys.argv[3:]
        stmts = []
        for source in sources:
            stmts += ac.load_statements(source)
        pkldump(stmts, output)
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
    elif cmd == 'filter_source':
        input_file = sys.argv[2]
        source_name = sys.argv[3]
        output_file = sys.argv[4]
        stmts = pklload(input_file)
        filt_stmts = [s for s in stmts
                      if s.evidence[0].source_api == source_name]
        ac.dump_statements(filt_stmts, prefixed_pkl(output_file))
    elif cmd == 'map_sequence':
        stmts = pklload(input_file)
        stmts = ac.map_sequence(stmts, save=prefixed_pkl(output_file))
    elif cmd == 'preassembled':
        stmts = pklload(input_file)
        cur_scorer = CuratedScorer()
        stmts = ac.run_preassembly(stmts, return_toplevel=False,
                                   belief_scorer=cur_scorer,
                                   save=prefixed_pkl(output_file),
                                   poolsize=16)

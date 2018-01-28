from util import *
import indra.tools.assemble_corpus as ac

if __name__ == '__main__':
    #sources = ['bel', 'biopax', 'reach', 'sparser']
    sources = [
        'output/bioexp_bel.pkl',
        'data/bioexp_biopax.pkl',
        'data/bioexp_sparser.pkl']
    stmts = []
    for source in sources:
        stmts += ac.load_statements(source)
    pkldump(stmts, 'all_raw')
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

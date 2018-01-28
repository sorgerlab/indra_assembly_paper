from util import pklload, listify_dict
import pickle
import indra.tools.assemble_corpus as ac

with open('korkut_genes.txt', 'rt') as fh:
    korkut_genes = [l.strip() for l in fh.readlines()]

with open('pmids_for_gene.pkl', 'rb') as fh:
    pmid_dict = pickle.load(fh)

korkut_pmids = {gene: pmid_dict[gene] for gene in korkut_genes}
korkut_pmids_list = listify_dict(korkut_pmids)

full_stmtsd = {}
full_stmts = {}
korkut_stmtsd = {}
korkut_stmts = {}
relevant_korkut_stmts = {}
relevant_full_stmts = {}
for reader in ['sparser', 'reach']:
    print('Running for %s' % reader)
    full_stmtsd[reader] = pklload(reader)
    korkut_stmtsd[reader] = {pmid: full_stmtsd[reader].get(pmid, [])
                             for pmid in korkut_pmids_list}

    full_stmts[reader] = listify_dict(full_stmtsd[reader]) 
    korkut_stmts[reader] = listify_dict(korkut_stmtsd[reader]) 

    full_stmts[reader] = ac.map_grounding(full_stmts[reader])
    korkut_stmts[reader] = ac.map_grounding(korkut_stmts[reader])

    relevant_korkut_stmts[reader] = \
        ac.filter_gene_list(korkut_stmts[reader],
                            korkut_genes, 'one')

    relevant_full_stmts[reader] = ac.filter_gene_list(
        full_stmts[reader], korkut_genes, 'one')

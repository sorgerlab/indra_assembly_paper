import pickle
import indra.tools.assemble_corpus as ac

with open('korkut_genes.txt', 'rt') as fh:
    korkut_genes = [l.strip() for l in fh.readlines()]

with open('pmids_for_gene.pkl', 'rb') as fh:
    pmid_dict = pickle.load(fh)

korkut_pmids = {gene: pmid_dict[gene] for gene in korkut_genes}
korkut_pmids_list = listify_dict(korkut_pmids)

stmts = {}
stmts_list = {}
korkut_stmts = {}
relevant_stmts_list = {}
for reader in ['sparser', 'reach']:
    with open('/home/bmg16/data/%s.pkl' % reader, 'rb') as fh:
        print('Running for %s' % reader)
        stmts[reader] = pickle.load(fh)
        korkut_stmts[reader] = {pmid: stmts[reader].get(pmid, [])
                                 for pmid in korkut_pmids_list}
        stmts_list[reader] = listify_dict(stmts[reader])
        stmts_list[reader] = ac.map_grounding(stmts_list[reader])
        relevant_stmts_list[reader] = ac.filter_gene_list(
            stmts_list[reader], korkut_genes, 'one')

from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import sys
import pickle
import logging
from random import shuffle
from indra.statements import Complex
from indra.preassembler import Preassembler
from indra.preassembler.hierarchy_manager import hierarchies
from indra.databases import biogrid_client as bg
from indra.databases import hgnc_client, uniprot_client
from indra.util import write_unicode_csv
from indra.tools import assemble_corpus as ac

logger = logging.getLogger('complexes')


def analyze(filename, sample_size=None):
    all_stmts = ac.load_statements(filename)

    complexes = ac.filter_by_type(all_stmts, Complex)
    protein_complexes = ac.filter_genes_only(complexes)
    protein_complexes = ac.filter_human_only(protein_complexes)

    # Optionally take a random sample
    if sample_size:
        shuffle(protein_complexes)
        protein_complexes = protein_complexes[0:sample_size]

    rows = [('Raw complexes', len(complexes)),
            ('Complexes between genes', len(protein_complexes))]

    logger.info('Mapping gene IDs to gene symbols')
    gene_ids = list(set([ag.db_refs['HGNC'] for stmt in protein_complexes
                                            for ag in stmt.members]))
    genes = [hgnc_client.get_hgnc_name(id) for id in gene_ids]

    # Get complexes from BioGrid and combine duplicates
    num_genes_per_query = 50
    start_indices = range(0, len(genes), num_genes_per_query)
    end_indices = [i + num_genes_per_query
                   if i + num_genes_per_query < len(genes) else len(genes)
                   for i in start_indices]
    bg_complexes = []
    for i in range(len(start_indices)):
        logger.info("Querying biogrid for %s" %
                    str(genes[start_indices[i]:end_indices[i]]))
        bg_complexes += (bg.get_statements(
                                genes[start_indices[i]:end_indices[i]]))

    # Filter out Biogrid statements not involving genes in the gene list
    # (this will make duplicate removal more efficient
    bg_filt = []
    for stmt in bg_complexes:
        if stmt.members[0].name in genes and \
           stmt.members[1].name in genes:
            bg_filt.append(stmt)
    # Might as well free up some memory
    del bg_complexes

    logger.info("Combining duplicates with biogrid...")
    pa = Preassembler(hierarchies, bg_filt + protein_complexes)
    pa.combine_duplicates()

    indra_only = []
    bg_only = []
    indra_and_bg = []
    for stmt in pa.unique_stmts:
        evidence_source_list = set([])
        for e in stmt.evidence:
            evidence_source_list.add(e.source_api)
        if evidence_source_list == set(['biogrid']):
            bg_only.append(stmt)
        elif not 'biogrid' in evidence_source_list:
            indra_only.append(stmt)
        else:
            indra_and_bg.append(stmt)

    stmt_header = ['STATEMENT_INDEX', 'AGENT_A', 'AGENT_B', 'NUM_MENTIONS',
                   'IN_BIOGRID']
    stmt_rows = [stmt_header]
    mention_header = ['STATEMENT_INDEX', 'AGENT_A', 'AGENT_B', 'READER',
                      'PMID', 'TEXT']
    mention_rows = [mention_header]
    index = 1
    for stmt in pa.unique_stmts:
        if stmt in bg_only:
            continue
        non_bg_ev_counts = len([e for e in stmt.evidence
                                if e.source_api != 'biogrid'])
        stmt_row = [index, stmt.members[0].name, stmt.members[1].name,
                    non_bg_ev_counts]
        mention_row = [index, stmt.members[0].name, stmt.members[1].name]
        if stmt in indra_and_bg:
            stmt_row.append(1)
            stmt_rows.append(stmt_row)
        elif stmt in indra_only:
            stmt_row.append(0)
            stmt_rows.append(stmt_row)
            for e in stmt.evidence:
                mention_rows.append(
                        mention_row + [e.source_api, e.pmid, e.text])
        else:
            assert False
        index += 1
    #for stmt in indra_only:
    #    rows.append([stmt.members[0].name, stmt.members[1].name,
    #                 str(len(stmt.evidence))])
    write_unicode_csv('complex_stmts.tsv', stmt_rows, delimiter='\t')
    write_unicode_csv('complex_stmt_mentions.tsv', mention_rows,
                      delimiter='\t')

    return {'indra_only': indra_only,
            'bg_only': bg_only,
            'indra_and_bg': indra_and_bg}


if __name__ == '__main__':
    # Load the statements
    if len(sys.argv) < 2:
        print("Usage: %s stmts_file" % sys.argv[0])
        sys.exit()
    results = analyze(sys.argv[1], sample_size=1000)
    print(results)


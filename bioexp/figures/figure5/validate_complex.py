from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import sys
import pickle
import logging
import random
from indra.statements import Complex
from indra.preassembler import Preassembler
from indra.preassembler.hierarchy_manager import hierarchies
from indra.databases import biogrid_client as bg
from indra.databases import hgnc_client, uniprot_client
from indra.util import write_unicode_csv
from indra.tools import assemble_corpus as ac

logger = logging.getLogger('complexes')

def preprocess_stmts(filename, sample_size=None):
    all_stmts = ac.load_statements(filename)

    complexes = ac.filter_by_type(all_stmts, Complex)
    # Require HGNC grounding
    hgnc_stmts = []
    for stmt in complexes:
        if all([('HGNC' in ag.db_refs) for ag in stmt.agent_list()]):
            hgnc_stmts.append(stmt)
    # Optionally take a random sample
    if sample_size:
        random.shuffle(hgnc_stmts)
        hgnc_stmts = hgnc_stmts[0:sample_size]
    return hgnc_stmts


def get_biogrid_stmts(stmts):
    logger.info('Mapping gene IDs to gene symbols')
    gene_ids = list(set([ag.db_refs['HGNC'] for stmt in stmts
                                            for ag in stmt.members]))
    genes = [hgnc_client.get_hgnc_name(id) for id in gene_ids]

    # Get complexes from BioGrid
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
    return bg_filt


def analyze_raw(filename, output_file, sample_size=None):
    protein_complexes = preprocess_stmts(filename, sample_size)
    bg_filt = get_biogrid_stmts(protein_complexes)
    # Build set of matches keys for Biogrid statements
    bg_keys = set([s.matches_key() for s in bg_filt])
    # Check each statement to see if there is a corresponding matches key
    # in the Biogrid set
    stmt_header = ['STATEMENT_INDEX',
                   'AGENT_A_NAME', 'AGENT_A_TEXT', 'AGENT_A_GROUNDING',
                   'AGENT_B_NAME', 'AGENT_B_TEXT', 'AGENT_B_GROUNDING',
                   'IN_BIOGRID', 'PMID', 'TEXT']
    rows = [stmt_header]

    def format_agent_entries(agent):
        db_refs_str = ','.join(['%s|%s' % (k, v)
                                for k, v in agent.db_refs.items()])
        ag_text = agent.db_refs.get('TEXT')
        text_str = ag_text if ag_text else ''
        return [agent.name, text_str, db_refs_str]

    for ix, stmt in enumerate(protein_complexes):
        ix += 1
        in_biogrid = 1 if stmt.matches_key() in bg_keys else 0
        if len(stmt.agent_list()) > 2:
            logger.info("Skipping statement with more than two members: %s"
                        % stmt)
            continue
        (ag_a, ag_b) = stmt.agent_list()
        row = [ix] + format_agent_entries(ag_a) + format_agent_entries(ag_b) +\
              [in_biogrid, stmt.evidence[0].pmid, stmt.evidence[0].text]
        rows.append(row)
    write_unicode_csv(output_file, rows, delimiter='\t')


def analyze_preassembled(filename, sample_size=None):
    protein_complexes = preprocess_stmts(filename, sample_size)
    bg_filt = get_biogrid_stmts(protein_complexes)

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
    write_unicode_csv('complex_stmts.tsv', stmt_rows, delimiter='\t')
    write_unicode_csv('complex_stmt_mentions.tsv', mention_rows,
                      delimiter='\t')

    return {'indra_only': indra_only,
            'bg_only': bg_only,
            'indra_and_bg': indra_and_bg}


if __name__ == '__main__':
    # Load the statements
    if len(sys.argv) < 3:
        print("Usage: %s stmts_file output_file" % sys.argv[0])
        sys.exit()
    random.seed(1)
    results = analyze_raw(sys.argv[1], sys.argv[2], sample_size=1000)

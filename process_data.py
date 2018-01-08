import pandas
from indra.databases import uniprot_client
from indra.databases import hgnc_client


data_file = 'data/Korkut et al. Data 05122017.xlsx'
ras227_file = 'data/ras_pathway_proteins.csv'
drug_grounding_file = 'data/drug_grounding.csv'


def read_data(fname=data_file):
    """Returns the data as a dictionary."""
    if fname.endswith('2017.xlsx'):
        skiprows1 = []
        skiprows2 = []
    else:
        skiprows1 = [0]
        skiprows2 = range(5)
    data = {}
    data['protein'] = pandas.read_excel(fname, sheetname='Protein Data',
                                        skiprows=skiprows1, index_col=None)
    data['phenotype'] = pandas.read_excel(fname,
                                          sheetname='Phenotype Data',
                                          skiprows=skiprows1, index_col=None)
    data['antibody'] = pandas.read_excel(fname,
                                         sheetname='Antibody Data',
                                         skiprows=skiprows2, index_col=None)
    data['prediction'] = pandas.read_excel(fname,
                                           sheetname='Prediction Targets',
                                           index_col=None)
    return data


def get_ras227_genes():
    df = pandas.read_csv(ras227_file, sep='\t', index_col=None, header=None,
                         encoding='utf-8')
    gene_names = [x for x in df[0]]
    return gene_names


def get_drug_targets(fname=None):
    if not fname:
        fname = drug_grounding_file
    df = pandas.read_csv(fname, index_col=None, header=None, encoding='utf-8')
    abbrevs = df[1]
    target_upids = df[6]
    targets = {}
    for abb, tupid in zip(abbrevs, target_upids):
        targets[abb] = [uniprot_client.get_gene_name(ui)
                        for ui in tupid.split(',')]
    return targets


def get_all_gene_names(data, out_file='prior_genes.txt'):
    """Return all gene names corresponding to all ABs."""
    filt = pandas.notnull(data['antibody']['Protein Data ID'])
    data_filt = data['antibody'][filt]
    gene_names = data_filt['Gene Name']
    uniprot_ids = data_filt['UniProt ID']
    all_genes = set()
    invalid_genes = set()
    for gn, upid in zip(gene_names, uniprot_ids):
        # Some entries are lists of genes separated by commas
        # and we also strip off extra spaces
        names = [x.strip() for x in gn.split(',')]
        ids = [x.strip() for x in upid.split(',')]
        names_from_ids = [uniprot_client.get_gene_name(x) for x in ids]
        # Find invalid gene names
        for name in names:
            if not hgnc_client.get_hgnc_id(name):
                print('Invalid or deprecated gene symbol: %s' % name)
                invalid_genes.add(name)
        # Find inconsistent gene names and UniProt IDs
        if set(names) != set(names_from_ids):
            print('Inconsistent entries:')
            print('- Given gene names: %s' % ','.join(names))
            print('- Genes from uniprot IDs: %s' % ','.join(names_from_ids))
        # Add both the gene names and the gene names derived from UniProt IDs
        all_genes = all_genes.union(set(names)).union(set(names_from_ids))
    # Finally remove the invalid gene names
    all_genes = list(all_genes.difference(invalid_genes))
    # Add drug target genes
    drug_targets = get_drug_targets()
    for targets in drug_targets.values():
        all_genes += targets
    # Add other important genes, for now, the RAS pathway
    all_genes += get_ras227_genes()
    all_genes = sorted(list(set(all_genes)))
    print('%d genes in total' % len(all_genes))
    with open(out_file, 'wb') as fh:
        for gene in all_genes:
            fh.write(('%s\n' % gene).encode('utf-8'))
    return all_genes

if __name__ == '__main__':
    data = read_data(data_file)
    gene_names = get_all_gene_names(data)


import sys
import pandas
from indra.databases import uniprot_client
from indra.databases import hgnc_client
from indra.statements import Agent, ModCondition



def read_data(fname):
    """Returns the data as a dictionary."""
    if fname.endswith('2017.xlsx'):
        skiprows1 = []
        skiprows2 = []
    else:
        skiprows1 = [0]
        skiprows2 = range(5)
    data = {}
    data['protein'] = pandas.read_excel(fname, sheet_name='Protein Data',
                                        skiprows=skiprows1, index_col=None)
    data['phenotype'] = pandas.read_excel(fname,
                                          sheet_name='Phenotype Data',
                                          skiprows=skiprows1, index_col=None)
    data['antibody'] = pandas.read_excel(fname,
                                         sheet_name='Antibody Data',
                                         skiprows=skiprows2, index_col=None)
    data['prediction'] = pandas.read_excel(fname,
                                           sheet_name='Prediction Targets',
                                           index_col=None)
    return data


def get_ras227_genes(ras227_file):
    df = pandas.read_csv(ras227_file, sep='\t', index_col=None, header=None,
                         encoding='utf-8')
    gene_names = [x for x in df[0]]
    return gene_names


def get_drug_targets(fname):
    df = pandas.read_csv(fname, index_col=None, header=None, encoding='utf-8')
    abbrevs = df[1]
    target_upids = df[6]
    targets = {}
    for abb, tupid in zip(abbrevs, target_upids):
        targets[abb] = [uniprot_client.get_gene_name(ui)
                        for ui in tupid.split(',')]
    return targets


def get_all_gene_names(data, ras_file, drug_file, out_file=None):
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
    drug_targets = get_drug_targets(drug_file)
    for targets in drug_targets.values():
        all_genes += targets
    # Add other important genes, for now, the RAS pathway
    all_genes += get_ras227_genes(ras_file)
    all_genes = sorted(list(set(all_genes)))
    print('%d genes in total' % len(all_genes))
    if out_file:
        with open(out_file, 'wb') as fh:
            for gene in all_genes:
                fh.write(('%s\n' % gene).encode('utf-8'))
    return all_genes


def get_phospho_antibody_map(fname):
    # First gather the annotations for the phosphosites
    df = pandas.read_csv(fname, index_col=None, sep=',', encoding='utf8')
    antibody_map = {}

    for _, row in df.iterrows():
        ps = row['phosphosite']
        sub_upid = row['SUB_ID']
        if not pandas.isnull(sub_upid):
            if sub_upid.find('-') != -1:
                sub_upid = sub_upid.split('-')[0]
            sub_hgnc_symbol = uniprot_client.get_gene_name(sub_upid)
            sub_hgnc = hgnc_client.get_hgnc_id(sub_hgnc_symbol)
        else:
            sub_hgnc_symbol = row['SUB_GENE']
            sub_hgnc_id = hgnc_client.get_hgnc_id(sub_hgnc_symbol)
            sub_upid = hgnc_client.get_uniprot_id(sub_hgnc_id)
            if sub_upid is None:
                continue
        sub = Agent(sub_hgnc_symbol,
                    db_refs={'UP': sub_upid,'HGNC': sub_hgnc})
        residue = row['Actual_site'][0]
        if len(row['Actual_site']) > 1:
            position = row['Actual_site'][1:]
        else:
            position = None
        mc = ModCondition('phosphorylation', residue, position)
        sub.mods = [mc]
        if ps in antibody_map:
            found = False
            for p in antibody_map[ps]:
                if p.name == sub.name and p.mods[0].residue == residue and \
                    p.mods[0].position == position:
                    found = True
                    break
            if not found:
                antibody_map[ps].append(sub)
        else:
            antibody_map[ps] = [sub]
    return antibody_map


def get_antibody_map(data, ab_file):
    phos_ab_map = get_phospho_antibody_map(ab_file)
    ab_map = {}
    for _, row in data['antibody'].iterrows():
        ab_name = row['Protein Data ID']
        if ab_name in phos_ab_map:
            continue
        upids = row['UniProt ID'].split(',')
        for upid in upids:
            hgnc_symbol = uniprot_client.get_gene_name(upid)
            hgnc_id = hgnc_client.get_hgnc_id(hgnc_symbol)
            target = Agent(hgnc_symbol,
                           db_refs={'UP': upid,'HGNC': hgnc_id})
            try:
                ab_map[ab_name].append(target)
            except KeyError:
                ab_map[ab_name] = [target]
    ab_map.update(phos_ab_map)
    return ab_map


if __name__ == '__main__':
    data_file, ras_file, drug_file, output_file = sys.argv[1:5]

    data = read_data(data_file)
    gene_names = get_all_gene_names(data, ras_file, drug_file,
                                    output_file)


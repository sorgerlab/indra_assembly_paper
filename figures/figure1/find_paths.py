import csv
import pickle
from os.path import dirname, abspath, join
import networkx as nx
from indra.util import read_unicode_csv

csv.field_size_limit(150000) # Accommodate a particularly long line

network_dir = join(dirname(abspath(__file__)), '..', '..', 'networks')


def load_pc_network(filter_genes=None, flatten=True):
    """Get Pathway Common network as a networkx MultiDiGraph.

    Parameters
    ----------
    filter_genes : list of str
        List of gene names to filter the network against. Only edges between
        nodes in this list will be added to the network.
    flatten : boolean
        Whether to flatten the MultiDiGraph (which permits multiple edges between
        the same nodes) to a DiGraph (single edge between two nodes).
        If flattened, edge metadata (e.g., interaction type) is discarded.
        Useful for running networkx pathfinding functions, which are not
        implemented for MultiDiGraphs. Default is True.
    """
    pc_filename = join(network_dir, 'PathwayCommons9.All.hgnc.txt')
    # Column names in the TSV file
    col_names = ['PARTICIPANT_A', 'INTERACTION_TYPE', 'PARTICIPANT_B',
                 'INTERACTION_DATA_SOURCE', 'INTERACTION_PUBMED_ID',
                 'PATHWAY_NAMES', 'MEDIATOR_IDS']
    # Get the data, skipping the header line
    print("Processing Pathway Commons TSV file")
    pc_data_generator = read_unicode_csv(pc_filename, delimiter='\t', skiprows=1)
    # Load into a networkx MultiDiGraph
    if flatten:
        pc_graph = nx.DiGraph()
    else:
        pc_graph = nx.MultiDiGraph()
    for ix, row in enumerate(pc_data_generator):
        # Handle possible missing rows
        if not row:
            continue
        if (ix+1) % 100000 == 0:
            print("Row %d" % (ix+1))
        subj = row[0]
        obj = row[2]
        # If desired, these lines can be uncommented to put more of the extended
        # SIF metadata into the network edges
        edge_data = dict(zip(col_names[3:-1], row[3:-1]))
        edge_data['relation'] = row[1]
        if not filter_genes or \
           (subj in filter_genes and obj in filter_genes):
            if flatten:
                pc_graph.add_edge(subj, obj)
            else:
                pc_graph.add_edge(subj, obj, attr_dict=edge_data)
    return pc_graph




if __name__ == '__main__':
    rerun = False
    if rerun:
        prior_genes_file = '../../prior_genes.txt'
        with open(prior_genes_file, 'rt') as f:
            prior_genes = [line.strip() for line in f.readlines()]
        pc_graph = load_pc_network(filter_genes=prior_genes)
        with open('pc_graph.pkl', 'wb') as f:
            pickle.dump(pc_graph, f)
    else:
        with open('pc_graph.pkl', 'rb') as f:
            pc_graph = pickle.load(f)

    source = 'BRAF'
    target = 'MAPK1'
    for path in nx.shortest_simple_paths(pc_graph, source, target):
        if len(path) > 3:
            break
        print(path)
    """
    canonical_path = ['EGFR', 'GRB2', 'SOS1', 'KRAS', 'BRAF', 'MAP2K1', 'MAPK1']
    for i in range(len(canonical_path)):
        head = canonical_path[i]
        for node in canonical_path[i+1:]:
            if pc_graph.has_edge(head, node):
                print("%s, %s" % (head, node))
    """


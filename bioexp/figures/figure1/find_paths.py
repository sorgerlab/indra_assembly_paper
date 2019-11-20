import sys
import csv
import pickle
from os.path import dirname, abspath, join
import networkx as nx
from indra.util import read_unicode_csv
import paths_graph as pg

csv.field_size_limit(150000) # Accommodate a particularly long line



def filter_direct(mdg):
    """Split MultiDiGraph into two DiGraphs, direct and indirect.

    Edges are considered direct if at least one of the edges has the
    relation 'interacts-with'.
    """
    direct_edges = set()
    all_edges = set()
    for u, v, edge_data in pc_graph.edges_iter(data=True):
        all_edges.add((u, v))
        if edge_data['relation'] in ('interacts-with', 'in-complex-with'):
            direct_edges.add((u, v))
    indirect_edges = all_edges - direct_edges
    pc_graph_dir = nx.DiGraph()
    pc_graph_indir = nx.DiGraph()
    pc_graph_dir.add_edges_from([(u, v, {'direct': True})
                                 for u, v in direct_edges])
    pc_graph_indir.add_edges_from([(u, v, {'direct': False})
                                   for u, v in indirect_edges])
    return (pc_graph_dir, pc_graph_indir)


def flatten_network(mdg):
    """Combined edges to convert MultiDiGraph to DiGraph."""
    all_edges = set((u, v) for u, v in pc_graph.edges_iter())
    dg = nx.DiGraph()
    dg.add_edges_from(all_edges)
    return dg


def load_pc_network(pc_filename, filter_genes=None):
    """Get Pathway Common network as a networkx MultiDiGraph.

    Parameters
    ----------
    filter_genes : list of str
        List of gene names to filter the network against. Only edges between
        nodes in this list will be added to the network.
    flatten : boolean
        Whether to flatten the MultiDiGraph (which permits multiple edges
        between the same nodes) to a DiGraph (single edge between two nodes).
        If flattened, edge metadata (e.g., interaction type) is discarded.
        Useful for running networkx pathfinding functions, which are not
        implemented for MultiDiGraphs. Default is True.
    """
    # Column names in the TSV file
    col_names = ['PARTICIPANT_A', 'INTERACTION_TYPE', 'PARTICIPANT_B',
                 'INTERACTION_DATA_SOURCE', 'INTERACTION_PUBMED_ID',
                 'PATHWAY_NAMES', 'MEDIATOR_IDS']
    # Get the data, skipping the header line
    print("Processing Pathway Commons TSV file")
    pc_data_generator = read_unicode_csv(pc_filename, delimiter='\t',
                                         skiprows=1)
    # Load into a networkx MultiDiGraph
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
            pc_graph.add_edge(subj, obj, attr_dict=edge_data)
    return pc_graph


if __name__ == '__main__':

    # Script to parse the PC data and cache the network
    if sys.argv[1] == 'parse_pc':
        pc_file = sys.argv[2]
        prior_genes_file = sys.argv[3]
        pc_pickle = sys.argv[4]
        with open(prior_genes_file, 'rt') as f:
            prior_genes = [line.strip() for line in f.readlines()]
        pc_graph = load_pc_network(pc_file)
        with open(pc_pickle, 'wb') as f:
            pickle.dump(pc_graph, f)
    # Script to work off of the cached pickle file and generate paths
    elif sys.argv[1] == 'find_paths':
        pc_pickle = sys.argv[2]
        source = 'RAF1'
        target = 'MAPK1'
        print("Loading pickle")
        with open(pc_pickle, 'rb') as f:
            pc_graph = pickle.load(f)
        pc_graph = flatten_network(pc_graph)
        max_depth = 6
        print("Counting paths...")
        import time
        start = time.time()
        f_reach, b_reach = pg.get_reachable_sets(pc_graph, source, target,
                                                 max_depth)
        for length in range(1, max_depth+1):
            print("Getting path of length %d between %s and %s" %
                  (length, source, target))
            pg_start = time.time()
            paths_graph = pg.PathsGraph.from_graph(pc_graph, source, target,
                                                   length, f_reach, b_reach)
            print("PathsGraph has %d nodes" % len(paths_graph.graph))
            pg_elapsed = time.time() - pg_start
            print("PG took time", pg_elapsed)
            cfpg_start = time.time()
            cfpg = pg.CFPG.from_pg(paths_graph)
            print("CFPG has %d nodes" % len(cfpg.graph))
            cfpg_elapsed = time.time() - cfpg_start
            print("CFPG took time", cfpg_elapsed)
            path_count = cfpg.count_paths()
            print("%d paths of length %d" % (path_count, length))
        elapsed = time.time() - start
        print("PG", elapsed)
        # Calculating time for NX
        start = time.time()
        paths = list(nx.all_simple_paths(pc_graph, source, target, max_depth))
        elapsed = time.time() - start
        print("NX", elapsed)
"""

        from collections import defaultdict
        path_counts = defaultdict(lambda: 0)
        last_length = 0
        for ix, path in enumerate(nx.shortest_simple_paths(pc_graph,
                                                           source, target)):
            path_len = len(path) - 1
            if (ix+1) % 1000 == 0:
                print("Path # %d, len %d" % ((ix+1), path_len))
                print(path_counts)
            if path_len > max_depth:
                break
            path_counts[path_len] += 1

        #pc_graph_dir, pc_graph_indir = filter_direct(pc_graph)
        paths = []
        for g in (pc_graph_dir,):
            for path in nx.shortest_simple_paths(g, source, target):
                if len(path) > 5:
                    break
                paths.append(path)
        paths_output_file = join(build_dir, 'fig1_pc_egfr_mapk1_paths.txt')
        with open(paths_output_file, 'wt') as f:
            for path in paths:
                f.write('%s\n' % str(path))

    canonical_path = ['EGFR', 'GRB2', 'SOS1', 'KRAS', 'BRAF', 'MAP2K1', 'MAPK1']
    for i in range(len(canonical_path)):
        head = canonical_path[i]
        for node in canonical_path[i+1:]:
            if pc_graph.has_edge(head, node):
                print("%s, %s" % (head, node))
"""


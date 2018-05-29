import sys
import random
from collections import defaultdict
import numpy as np
from matplotlib import pyplot as plt
from matplotlib_venn import venn3
from indra.tools import assemble_corpus as ac
from util import pklload


def plot_statement_overlap(stmts, plot_filename):
    """Generate a Venn diagram showing reader overlap (for REACH, Medscan, and
    Sparesr) among duplicate statements."""
    # Iterate over the preassembled statements, collecting lists of UUIDs
    # with evidence from each source
    sources = defaultdict(set)
    for stmt in stmts:
        for ev in stmt.evidence:
            sources[ev.source_api].add(stmt.uuid)
    # Plot the Venn diagram
    plt.figure()
    subsets = (sources['reach'], sources['medscan'], sources['sparser'])
    venn3(subsets=subsets, set_labels=('REACH', 'Medscan', 'Sparser'))
    plt.title('Statement overlap among readers')
    plt.savefig(plot_filename)


def plot_belief_distributions(stmts_dict, basename):
    """For each reader and for the combined reading only results, make a bar
    plot of the distribution of belief scores in different bins. In addition,
    calculate the expected number of statements in each bin if there were
    not overlap among readers, and compare to the observed numbers from the
    combined results from all readers.
    """
    bins = [0.0, 0.5, 0.8, 0.9, 0.99, 1.0]
    def plot_scores(source, belief_scores):
        hist_result = plt.hist(scores, bins)
        counts = hist_result[0]
        plt.figure()
        plt.bar(range(len(counts)), counts,
                tick_label=('<0.5', '0.5-0.8', '0.8-0.9', '0.9-0.99', '>0.99'))
        plt.xlabel('Belief score')
        plt.ylabel('Number of Statements')
        ax = plt.gca()
        ax.set_xticklabels(('< 0.5', '0.5-0.8', '0.8-0.9', '0.9-0.99',
                            '> 0.99'))
        plt.savefig('%s_%s.pdf' % (basename, source))
        return counts

    counts_by_source = {}
    totals = np.zeros(len(bins) - 1)
    for source, stmts in stmts_dict.items():
        scores = [s.belief for s in stmts]
        counts = plot_scores(source, scores)
        counts_by_source[source] = counts
        if source != 'reading':
            totals += counts
    # Print counts for each bin for each source
    print("Belief score counts by source")
    print("Bins: %s" % str(bins))
    print(counts_by_source)
    print("Totals for each bin across each source")
    print(totals)

    # Now, make a final plot comparing totals vs. assembled
    def plot_belief_comparison(values1, values2, type):
        if type not in ('percent', 'count'):
            raise ValueError('type must be one of "percent", "count"')
        bins = [0.0, 0.5, 0.8, 0.9, 0.99, 1.0]
        plt.figure()
        index = np.array(range(len(total_counts)))
        width = 0.35
        plt.bar(index, values1, width=width)
        plt.bar(index+width, values2, width=width)
        plt.xlabel('Belief score')
        if type == 'percent':
            plt.ylabel('Pct. of Statements')
        elif type == 'count':
            plt.ylabel('Number of Statements')
        ax = plt.gca()
        ax.set_xticks(index + width/2)
        ax.set_xticklabels(('< 0.5', '0.5-0.8', '0.8-0.9', '0.9-0.99',
                            '> 0.99'))
        plt.subplots_adjust(left=0.16)
        plt.show()
        plt.savefig('%s_comparison_%s.pdf' % (basename, type))

    read_pcts = (100 * counts_by_source['reading'] /
                 np.sum(counts_by_source['reading'])
    total_pcts = 100 * totals / np.sum(totals)

    plot_belief_comparison(totals, counts_by_source['reading'], 'count')
    plot_belief_comparison(total_pcts, read_pcts, 'percent')


if __name__ == '__main__':
    # Load statements
    stmts_dict = {}
    stmts_dict['reach'] = pklload('reach_only_preassembled')
    stmts_dict['medscan'] = pklload('medscan_only_preassembled')
    stmts_dict['sparser'] = pklload('sparser_only_preassembled')
    stmts_dict['reading'] = pklload('reading_only_preassembled')

    # Make plots
    plot_statement_overlap(stmts_dict['reading'],
                           'output/stmt_overlap_reading.pdf')
    plot_belief_distributions(stmts_dict, 'output/belief_scores')

    # Sample statements from each belief bin
    stmt_bins = ((0.5, 0.8, '0.5_0.8'),
                 (0.8, 0.9, '0.8_0.9'),
                 (0.9, 0.99, '0.9_0.99'),
                 (0.99, 1.0, '0.9_1.0'),
                )
    random.seed(1)
    for lbound, ubound, label in stmt_bins:
        stmts_by_belief = [s for s in stmts_dict['reading']]
        random.shuffle(stmts_by_belief)
        pkldump(stmts_by_belief[0:1000], 'reading_belief_%s_sample' % label)


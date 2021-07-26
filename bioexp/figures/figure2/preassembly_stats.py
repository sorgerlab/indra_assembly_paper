import pickle
from collections import Counter
from os.path import dirname, join
from matplotlib import pyplot as plt
from indra.preassembler import render_stmt_graph
from bioexp.util import set_fig_params, fontsize, format_axis, red, based

def plot_frequencies(counts, x_label, y_label, fig_filename, plot_type='dot',
                     log_x=False, log_y=False):
    if plot_type not in ('dot', 'bar'):
        raise ValueError("plot_type must be one of ('bar', 'dot')")
    ctr = Counter(counts)
    ctr = sorted([(k, v) for k, v in ctr.items()],
                 key=lambda x: x[1], reverse=True)
    counts, stmts_per_count = zip(*ctr)
    fig_path = join(build_dir, fig_filename)
    fig = plt.figure(figsize=(1.5, 1.5), dpi=150)
    if plot_type == 'dot':
        plt.plot(counts, stmts_per_count, linestyle='', marker='.',
                 markersize='1', color=red)
    elif plot_type == 'bar':
        plt.bar(counts, stmts_per_count, color=red)
    ax = fig.gca()
    if log_x:
        ax.set_xscale('log')
    if log_y:
        ax.set_yscale('log')
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    format_axis(ax, tick_padding=2)
    plt.subplots_adjust(left=0.22, bottom=0.18, right=0.94, top=0.97)
    plt.savefig(fig_path)


def depth_of_support(stmt, depth):
    if not stmt.supported_by:
        return 0
    else:
        depths = [depth_of_support(supp_stmt, depth)
                  for supp_stmt in stmt.supported_by]
        return max(depths) + 1


def render_stmt_support(stmts):
    for ix, stmt in enumerate(stmts):
        g = render_stmt_graph([stmt], english=True, reduce=True)
        fig_filename = 'fig2_stmt_graph_%d.pdf' % (ix+1)
        fig_path = join(build_dir, fig_filename)
        g.draw(fig_path, prog='dot')


if __name__ == '__main__':
    data_dir = join(dirname(__file__), '..', '..', '..', 'data')
    build_dir = based

    stmts_file = join(data_dir, 'bioexp_asmb_preassembled.pkl')

    # Load the pickle
    print("Loading statements from %s" % stmts_file)
    with open(stmts_file, 'rb') as f:
        stmts = pickle.load(f)
    print("%d stmts" % len(stmts))

    set_fig_params()

    # Get/plot evidence distribution
    ev_counts = [len(s.evidence) for s in stmts]
    plot_frequencies(ev_counts, 'Number of evidences',
                     'Number of statements', 'fig2_evidence_distribution.pdf',
                     plot_type='dot', log_x=True, log_y=True)

    # Supported-by distribution
    supp_by_counts = [len(s.supported_by) for s in stmts]
    plot_frequencies(supp_by_counts, 'Supporting statements',
                     'Number of statements',
                     'fig2_supported_by_distribution.pdf',
                     plot_type='bar', log_y=True)

    # Get depths of support for each statement
    supp_depths_stmts = [(s, depth_of_support(s, 0)) for s in stmts]
    supp_depths_stmts = sorted(supp_depths_stmts, key=lambda x: x[1],
                               reverse=True)
    supp_depths = [t[1] for t in supp_depths_stmts]
    plot_frequencies(supp_depths, 'Depth of support', 'Number of statements',
                     'fig2_depths_of_support.pdf',
                     plot_type='bar', log_y=True)
    num_graphs = 30
    render_stmt_support([t[0] for t in supp_depths_stmts[0:num_graphs]])

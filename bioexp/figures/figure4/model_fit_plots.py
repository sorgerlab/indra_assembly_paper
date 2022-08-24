import sys
import pickle
from os.path import join
from matplotlib import pyplot as plt
from bioexp.util import set_fig_params, format_axis, fontsize

if __name__ == '__main__':
    fit_results_path = sys.argv[1]
    reader = sys.argv[2]
    output_dir = sys.argv[3]

    with open(fit_results_path, 'rb') as f:
        results = pickle.load(f)

    # Get results entries
    bin_res = results[3]
    betabin_res = results[5]
    bel_res = results[1]

    # Plot binomial
    ax = bin_res[1].plot_stmt_fit(bin_res[2], 'Binomial', color='r')
    # Plot beta-binomial
    ax = betabin_res[1].plot_stmt_fit(betabin_res[2], 'Beta-binomial',
                                      color='g', ax=ax)
    # Plot belief
    ax = bel_res[1].plot_stmt_fit(bel_res[2], 'Belief', color='b', ax=ax)
    plt.legend(frameon=False, loc='lower right', fontsize=fontsize)

    format_axis(ax)
    plt.subplots_adjust(left=0.152)
    plt.savefig(join(output_dir, f'fig4_{reader}_model_fits.pdf'))

    # Plot evidence distribution fits
    betabin_res[1].plot_ev_fit(betabin_res[2],
                           f'{reader.upper()} statements, Beta-binomial model')
    plt.savefig(join(output_dir, f'fig4_{reader}_ev_fit_betabin.pdf'))

    # Plot evidence distribution fits
    bel_res[1].plot_ev_fit(bel_res[2],
                           f'{reader.upper()} statements, Belief model')
    plt.savefig(join(output_dir, f'fig4_{reader}_ev_fit_belief.pdf'))

    # Plot evidence distribution fits
    bin_res[1].plot_ev_fit(bin_res[2],
                           f'{reader.upper()} statements, Binomial model')
    plt.savefig(join(output_dir, f'fig4_{reader}_ev_fit_bin.pdf'))


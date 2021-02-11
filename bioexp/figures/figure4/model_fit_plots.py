import pickle
from os.path import join
from matplotlib import pyplot as plt
from bioexp.util import set_fig_params, format_axis, fontsize

if __name__ == '__main__':
    fit_results_path = sys.argv[1]
    output_dir = sys.argv[2]

    with open(fit_results_path, 'rb') as f:
        results = pickle.load(f)

    # Get results entries
    bin_res = results[3]
    betabin_res = results[5]
    bel_res = results[1]

    plt.ion()

    # Plot binomial
    ax = bin_res[1].plot_stmt_fit(bin_res[2], 'Binomial', color='r')
    # Plot beta-binomial
    ax = betabin_res[1].plot_stmt_fit(betabin_res[2], 'Beta-binomial', color='g',
                                      ax=ax)
    # Plot belief
    ax = bel_res[1].plot_stmt_fit(bel_res[2], 'Belief', color='b', ax=ax)
    plt.legend(frameon=False, loc='lower right', fontsize=fontsize)

    format_axis(ax)
    plt.savefig(join(output_dir, 'fig4_reach_model_fits.pdf'))

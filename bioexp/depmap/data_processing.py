from os import environ
from os.path import join

import pystow
import numpy as np
import pandas as pd
from scipy import stats

from depmap_analysis.util.statistics import get_z, get_logp, get_n

depmap_release = '21q2'

# Base folder is either the BASEPATH environment variable or
# a folder under the pystow home, e.g., ~/.data/depmap_analysis/depmap/21q2
basedir_default = pystow.join('depmap_analysis', 'depmap', depmap_release).as_posix()
basedir = environ.get('BASEPATH', basedir_default)

### INPUT FILES
# http://www.broadinstitute.org/ftp/distribution/metabolic/papers/Pagliarini/MitoCarta3.0/Human.MitoCarta3.0.xls
mitocarta_file = join(basedir, 'Human.MitoCarta3.0.xls')

# https://ndownloader.figshare.com/files/27902376
sample_info_file = join(basedir, 'sample_info.csv')

# Download from https://ndownloader.figshare.com/files/13515395
rnai_file = join(basedir, 'D2_combined_gene_dep_scores.csv')

# Download from https://ndownloader.figshare.com/files/27902226
crispr_file = join(basedir, 'CRISPR_gene_effect.csv')

####


def get_corr(recalculate, data_df, filename):
    filepath = join(basedir, filename)
    if recalculate:
        data_corr = data_df.corr()
        data_corr.to_hdf('%s.h5' % filepath, filename)
    else:
        data_corr = pd.read_hdf('%s.h5' % filepath)
    return data_corr


def unstack_corrs(df):
    df_ut = df.where(np.triu(np.ones(df.shape), k=1).astype(np.bool))
    stacked: pd.DataFrame = df_ut.stack(dropna=True)
    return stacked


def filt_mitocorrs(ser):
    filt_ix = []
    filt_vals = []
    for (ix0, ix1), logp in ser.iteritems():
        if ix0 in mitogenes and ix1 in mitogenes:
            continue
        filt_ix.append((ix0, ix1))
        filt_vals.append(logp)
    index = pd.MultiIndex.from_tuples(filt_ix, names=['geneA', 'geneB'])
    filt_ser = pd.Series(filt_vals, index=index)
    return filt_ser


if __name__ == '__main__':
    recalculate = True

    # Process Mitocarta data
    mitocarta = pd.read_excel(mitocarta_file, sheet_name=1)
    mitogenes = list(mitocarta.Symbol.values)

    # Process cell line info from DepMap
    cell_line_df = pd.read_csv(sample_info_file)
    cell_line_map = cell_line_df[['DepMap_ID', 'CCLE_Name']]
    cell_line_map.set_index('CCLE_Name', inplace=True)

    # Process RNAi data
    rnai_df = pd.read_csv(rnai_file, index_col=0)
    rnai_df = rnai_df.transpose()
    gene_cols = ['%s' % col.split(' ')[0] for col in rnai_df.columns]
    rnai_df.columns = gene_cols
    rnai_df = rnai_df.join(cell_line_map)
    rnai_df = rnai_df.set_index('DepMap_ID')
    # Drop duplicate columns
    rnai_df = rnai_df.loc[:, ~rnai_df.columns.duplicated()]

    rnai_corr = get_corr(recalculate, rnai_df, 'rnai_correlations')
    rnai_n = get_n(recalculate, rnai_df, join(basedir, 'rnai_n'))
    rnai_logp = get_logp(recalculate, rnai_n, rnai_corr, join(basedir, 'rnai_logp'))
    rnai_z = get_z(recalculate, rnai_logp, rnai_corr, join(basedir, 'rnai_z_log'))

    # Process CRISPR data
    crispr_df = pd.read_csv(crispr_file, index_col=0)
    gene_cols = ['%s' % col.split(' ')[0] for col in crispr_df.columns]
    crispr_df.columns = gene_cols
    # Drop any duplicate columns (shouldn't be any for CRISPR, but just in case)
    crispr_df = crispr_df.loc[:, ~crispr_df.columns.duplicated()]

    crispr_corr = get_corr(recalculate, crispr_df, 'crispr_correlations')
    crispr_n = get_n(recalculate, crispr_df, join(basedir, 'crispr_n'))
    crispr_logp = get_logp(recalculate, crispr_n, crispr_corr, join(basedir, 'crispr_logp'))
    crispr_z = get_z(recalculate, crispr_logp, crispr_corr, join(basedir, 'crispr_z_log'))

    # Combine z-scores
    dep_z = (crispr_z + rnai_z) / np.sqrt(2)
    dep_z = dep_z.dropna(axis=0, how='all').dropna(axis=1, how='all')

    dep_logp = pd.DataFrame(np.log(2) + stats.norm.logcdf(-dep_z.abs()),
                            index=dep_z.columns, columns=dep_z.columns)

    filename = 'dep_z'
    z_filepath = join(basedir, '%s.h5' % filename)
    dep_z.to_hdf(z_filepath, filename)

    filename = 'dep_logp'
    logp_filepath = join(basedir, '%s.h5' % filename)
    dep_logp.to_hdf(logp_filepath, filename)

    df_logp = dep_logp
    total_comps = np.triu(~df_logp.isna(), k=1).sum()
    mitogenes_in_df = set(df_logp.columns).intersection(set(mitogenes))
    mito_comps = len(mitogenes_in_df)**2
    num_comps = total_comps - mito_comps

    alpha = 0.05
    bc_thresh = np.log(alpha / num_comps)
    sig_no_corr = unstack_corrs(df_logp[df_logp < np.log(alpha)])
    filt_corrs = filt_mitocorrs(sig_no_corr)
    sig_sorted = filt_corrs.sort_values().to_frame('logp')
    sig_sorted['rank'] = sig_sorted.rank()
    sig_sorted['bc_cutoff'] = bc_thresh
    sig_sorted['bh_crit_val'] = np.log((sig_sorted['rank'] / num_comps) * alpha)
    cm = np.log(num_comps) + np.euler_gamma + (1 / (2 * num_comps))

    sig_sorted['by_crit_val'] = sig_sorted['bh_crit_val'] - np.log(cm)
    filename = join(basedir, 'dep_stouffer_signif.pkl')
    sig_sorted.to_pickle(filename)

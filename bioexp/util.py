import os
import json
import time
import pickle
import matplotlib
from os.path import dirname, abspath, join


# CREATE A JSON FILE WITH THIS INFORMATION, E.G., a file consisting of:
# {"basename": "fallahi_eval", "basedir": "output"}
config_file = join(dirname(abspath(__file__)), '..', 'config.json')
with open(config_file, 'rt') as f:
    config = json.load(f)
# This is the base name used for all files created/saved
basen = config['basename']
# This is the base folder to read/write (potentially large) files from/to
# MODIFY ACCORDING TO YOUR OWN SETUP
based = config['basedir']


fontsize=6
red = '#e41a1c'
blue = '#377eb8'
green = '#4daf4a'
purple = '#984ea3'


# This makes it easier to make standardized pickle file paths
def prefixed_pkl(suffix):
    """Return full path to a pickle file based on name suffix."""
    return prefixed_file(suffix, 'pkl')


def prefixed_file(suffix, extension):
    """Return full path to a file based on name suffix and extension."""
    return os.path.join(based, basen + '_' + suffix + '.' + extension)


def pkldump(content, suffix):
    """Dump content into a pickle file based on a name suffix"""
    fname = prefixed_pkl(suffix)
    with open(fname, 'wb') as fh:
        pickle.dump(content, fh)


def pklload(suffix):
    """Load a pickle file based on a name suffix"""
    fname = prefixed_pkl(suffix)
    print('Loading %s' % fname)
    ts = time.time()
    with open(fname, 'rb') as fh:
        content = pickle.load(fh)
    te = time.time()
    print('Loaded %s in %.1f seconds' % (fname, te-ts))
    return content


def listify_dict(d):
    """Return a list by taking the union of dict entries that are lists."""
    ll = []
    for v in d.values():
        ll += v
    return list(set(ll))


def set_fig_params():
    matplotlib.rcParams['font.sans-serif'] = 'Arial'
    matplotlib.rcParams['text.usetex'] = True
    matplotlib.rcParams['text.latex.preamble'] = [
            '\\usepackage{helvet}',
            '\\usepackage{sansmath}',
            '\\sansmath',
            '\\usepackage{underscore}',]

def format_axis(ax, label_padding=2, tick_padding=0, yticks_position='left'):
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position(yticks_position)
    ax.yaxis.set_tick_params(which='both', direction='out', labelsize=fontsize,
                             pad=tick_padding, length=2, width=0.5)
    ax.xaxis.set_tick_params(which='both', direction='out', labelsize=fontsize,
                             pad=tick_padding, length=2, width=0.5)
    ax.xaxis.labelpad = label_padding
    ax.yaxis.labelpad = label_padding
    ax.xaxis.label.set_size(fontsize)
    ax.yaxis.label.set_size(fontsize)


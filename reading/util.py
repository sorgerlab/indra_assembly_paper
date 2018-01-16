from indra.util import _require_python3
import os
import json
import time
import pickle

# CREATE A JSON FILE WITH THIS INFORMATION, E.G., a file consisting of:
# {"basename": "fallahi_eval", "basedir": "output"}
with open('config.json', 'rt') as f:
    config = json.load(f)
# This is the base name used for all files created/saved
basen = config['basename']
# This is the base folder to read/write (potentially large) files from/to
# MODIFY ACCORDING TO YOUR OWN SETUP
based = config['basedir']


# This makes it easier to make standardized pickle file paths
def prefixed_pkl(suffix):
    """Return full path to a pickle file based on name suffix."""
    return os.path.join(based, basen + '_' + suffix + '.pkl')


def pkldump(suffix, content):
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

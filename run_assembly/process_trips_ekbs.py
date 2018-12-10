"""This script processes a set of TRIPS EKBs into Statements and saves
them as a pickle file."""
import glob
import pickle
from indra.sources import trips

fnames = glob.glob('../data/bioexp/trips/*.ekb')
print('Found %d files' % len(fnames))
all_stmts = []
for fname in fnames:
    tp = trips.process_xml_file(fname)
    all_stmts += tp.statements
with open('bioexp_trips.pkl', 'wb') as fh:
    pickle.dump(all_stmts, fh)

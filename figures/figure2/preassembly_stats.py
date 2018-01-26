from os.path import dirname, join
import pickle

data_dir = join(dirname(__file__), '..', '..', 'data')
build_dir = join(dirname(__file__), '..', '..', 'build')

stmts_file = join(data_dir, 'bioexp_bel.pkl')


# Load the pickle
with open(stmts_file, 'rb') as f:
    stmts = pickle.load(f)

# Write the number of statements to a file
output_file = join(build_dir, 'fig2_num_statements.txt')
with open(output_file, 'wt') as f:
    f.write('%d' % len(stmts))

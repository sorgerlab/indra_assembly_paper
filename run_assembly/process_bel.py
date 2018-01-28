from util import pkldump
from indra.sources import bel

bel_corpus = '../data/large_corpus-20170611.bel'

bp = bel.process_belscript(bel_corpus)
pkldump(bp.statements, 'bel')

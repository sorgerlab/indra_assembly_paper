from indra.sources import biogrid
from util import pkldump

bp = biogrid.BiogridProcessor()
pkldump(bp.statements, 'biogrid')


"""This script fixes a minor format issue with Conversion
Statements from BioPAX in which the list of incoming and outgoing
chemicals in the conversion can sometimes have lists in them.
This script flattens these lists."""

from util import pklload, pkldump
from indra.sources.biopax.processor import _list_listify
import indra.tools.assemble_corpus as ac
from indra.statements import Conversion

def merge(list_of_lists):
    ll = []
    for l in list_of_lists:
        ll += l
    return ll

# Load the pickle from the data directory
stmts = ac.load_statements('data/bioexp_biopax.pkl')
new_stmts = []
for stmt in stmts:
    if isinstance(stmt, Conversion):
        if any([isinstance(l, list) for l in stmt.obj_from]):
            print(stmt)
        stmt.obj_from = merge(_list_listify(stmt.obj_from))
        if any([isinstance(l, list) for l in stmt.obj_from]):
            print(stmt)
        stmt.obj_to = merge(_list_listify(stmt.obj_to))
stmts = pkldump(stmts, 'biopax_fixed')

import sys
from indra.assemblers.indranet import IndraNetAssembler
from bioexp.util import pklload, pkldump

if __name__ == '__main__':
    if len(sys.argv) != 2:
        usage = f"Usage: {sys.argv[0]} <pkl_basename>"
        print(usage)
        sys.exit(1)

    pkl_basename = sys.argv[1]
    stmts =  pklload(pkl_basename)
    ina = IndraNetAssembler(stmts)
    inet = ina.make_model()
    pkldump(inet, f"{pkl_basename}_indranet")

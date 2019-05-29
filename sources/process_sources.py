import os
import sys
import gzip
import glob
import indra
import shutil
import pickle
import zipfile
import logging
import urllib.request
from indra.sources import reach
from bioexp.transfer_s3 import download_from_s3, upload_to_s3


logger = logging.getLogger(__name__)


def process_source(source, cached, data_folder, target_folder):
    logger.info('Processing %s' % source)
    key = 'bioexp_%s.pkl' % source
    if cached:
        logger.info('Loading %s from cache' % key)
        download_from_s3(key, target_folder)
    else:
        fun = globals()['process_%s' % source]
        stmts = fun(data_folder)
        logger.info('Dumping %d statements into %s' % (len(stmts), key))
        fname = os.path.join(target_folder, key)
        with open(fname, 'wb') as fh:
            pickle.dump(stmts, fh)
        upload_to_s3(fname)


def process_pathway_commons(data_folder):
    # Reetreive the OWL file and put it in the data folder
    url = 'https://www.pathwaycommons.org/archives/PC2/v11/' + \
        'PathwayCommons11.All.BIOPAX.owl.gz'
    fname = os.path.join(data_folder, 'PathwayCommons11.All.BIOPAX.owl')
    logger.info('Downloading %s and extracting into %s' % (url, fname))
    gz_file = os.path.join(data_folder, 'PathwayCommons11.All.BIOPAX.owl.gz')
    urllib.request.urlretrieve(url, gz_file)
    with gzip.open(gz_file, 'rb') as fin:
        with open(fname, 'wb') as fout:
            shutil.copyfileobj(fin, fout)

    # Now process the OWL file to Statements
    from indra.sources import biopax
    bp = biopax.process_owl(fname)

    # Now filter out phosphosite
    stmts = [s for s in bp.statements if
             (s.evidence[0].annotations.get('source_sub_id')
              == 'phosphositeplus')]
    return stmts


def process_bel(data_folder):
    from indra.sources import bel
    url = 'https://arty.scai.fraunhofer.de/artifactory/bel/knowledge/' + \
        'large_corpus/large_corpus-20170611.bel'
    fname = os.path.join(data_folder, 'large_corpus-20170611.bel')
    urllib.request.urlretrieve(url, fname)
    bp = bel.process_belscript(fname)
    return bp.statements


def process_signor(data_folder):
    from indra.sources import signor
    sp = signor.process_from_web()
    return sp.statements


def process_rlimsp(data_folder):
    from indra.sources import rlimsp
    rlimsp_url = 'https://hershey.dbi.udel.edu/textmining/export/'
    medline_file = 'rlims.medline.json'
    pmc_file = 'rlims.pmc.json'

    stmts = []
    for fname, id_type in ((medline_file, 'pmid'), (pmc_file, 'pmcid')):
        out_file = os.path.join(data_folder, fname)
        url = rlimsp_url + fname
        urllib.request.urlretrieve(url, out_file)
        rp = rlimsp.process_from_json_file(out_file, id_type)
        stmts += rp.statements
    return stmts


def process_biogrid(data_folder):
    from indra.sources import biogrid
    bp = biogrid.BiogridProcessor()
    return bp.statements


def process_phosphosite(data_folder):
    from indra.sources import biopax
    fname = os.path.join(data_folder, 'Kinase_substrates.owl')
    bp = biopax.process_owl(fname)
    return bp.statements


def process_cbn(data_folder):
    from indra.sources import bel
    url = 'http://causalbionet.com/Content/jgf_bulk_files/Human-2.0.zip'
    zip_file = os.path.join(data_folder, 'Human-2.0.zip')
    urllib.request.urlretrieve(url, zip_file)
    cbn_folder = os.path.join(data_folder, 'cbn')
    os.mkdir(cbn_folder)
    with zipfile.ZipFile(zip_file) as fh:
        fh.extractall(path=cbn_folder)
    stmts = []
    for fname in glob.glob(os.path.join(cbn_folder, '*.jgf')):
        bp = bel.process_cbn_jgif_file(fname)
        stmts += bp.statements
    return stmts


def process_trrust(data_folder):
    from indra.sources import trrust
    tp = trrust.process_from_web()
    return tp.statements


def _process_reach_pmid(pmid):
    from indra.literature import s3_client
    try:
        logger.info('Processing %s' % pmid)
        reach_json_str = s3_client.get_reader_json_str('reach', pmid)
        rp = reach.process_json_str(reach_json_str, citation=pmid)
        return rp.statements
    except Exception as e:
        return []


def process_reach(data_folder):
    from multiprocessing import Pool
    pmid_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             os.pardir, 'run_assembly', 'pmids.txt')
    pmids = [l.strip() for l in open(pmid_file).readlines()]
    pool = Pool(4)
    stmts_ll = pool.map(_process_reach_pmid, pmids)
    pool.close()
    pool.join()
    stmts = []
    for stmts_l in stmts_ll:
        stmts += stmts_l
    return stmts


if __name__ == '__main__':
    data_folder = sys.argv[1]
    target_folder = sys.argv[2]
    cached = True if sys.argv[3] == 'True' else False
    sources = sys.argv[4:]
    print(sys.argv)
    for source in sources:
        process_source(source, cached, data_folder, target_folder)
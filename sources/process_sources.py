import os
import sys
import gzip
import json
import glob
import indra
import shutil
import pickle
import zipfile
import tarfile
import logging
import urllib.request
from indra.sources import reach, trips, medscan, sparser
from bioexp.transfer_s3 import download_from_s3, upload_to_s3


logger = logging.getLogger(__name__)
pmid_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         os.pardir, 'run_assembly', 'pmids.txt')


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


def process_hprd(data_folder):
    from indra.sources import hprd
    tgz_fname = 'HPRD_FLAT_FILES_041310.tar.gz'
    url = 'http://www.hprd.org/RELEASE9/%s' % tgz_fname
    local_tgz_fname = os.path.join(data_folder, tgz_fname)
    urllib.request.urlretrieve(url, local_tgz_fname)
    with tarfile.open(local_tgz_fname, 'r:gz') as fh:
        fh.extractall(data_folder)
    hprd_path = os.path.join(data_folder, 'FLAT_FILES_072010')
    hp = hprd.process_flat_files(
        id_mappings_file=os.path.join(hprd_path,
                                      'HPRD_ID_MAPPINGS.txt'),
        complexes_file=os.path.join(hprd_path,
                                    'PROTEIN_COMPLEXES.txt'),
        ptm_file=os.path.join(hprd_path,
                              'POST_TRANSLATIONAL_MODIFICATIONS.txt'),
        ppi_file=os.path.join(hprd_path,
                              'BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt'),
        seq_file=os.path.join(hprd_path,
                              'PROTEIN_SEQUENCES.txt')
        )
    return hp.statements



def _process_medscan_get_stmts(fname):
    mp = medscan.process_file(fname)
    if mp:
        return mp.statements
    else:
        return []


def process_medscan(data_folder):
    # NOTE: this function has not been run as is, unexpected issues could
    # come up
    from multiprocessing import Pool
    # Process PMID file first
    pmid_file = os.path.join(data_folder, 'medscan', 'DARPAcorpus.csxml')
    mp = medscan.process_file(pmid_file)
    pmid_stmts = mp.statements
    # Process PMC files next
    file_pattern = os.path.join(data_folder, 'medscan', 'pmids', '*.csxml')
    fnames = glob.glob(file_pattern)
    # Now process files in parallel
    pool = Pool(4)
    stmts_ll = pool.map(_process_medscan_get_stmts, fnames)
    # Flatten the list of lists
    stmts = []
    for sts in stmts_ll:
        stmts += sts
    # Only add PMID statements if their PMID isn't covered by PMC
    pmc_pmids = {s.evidence[0].pmid for s in stmts}
    for stmt in pmid_stmts:
        if stmt.evidence[0].pmid not in pmc_pmids:
            stmts.append(stmt)
    return stmts


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
    pmids = [l.strip() for l in open(pmid_file).readlines()]
    pool = Pool(4)
    stmts_ll = pool.map(_process_reach_pmid, pmids)
    pool.close()
    pool.join()
    stmts = []
    for stmts_l in stmts_ll:
        stmts += stmts_l
    return stmts


def _process_trips_fname(fname):
    tp = trips.process_xml_file(fname)
    if tp:
        pmid = os.path.basename(fname)[:-4]
        for stmt in tp.statements:
            for ev in stmt.evidence:
                ev.pmid = pmid
        return tp.statements
    else:
        return []


def process_trips(data_folder):
    from multiprocessing import Pool
    file_pattern = os.path.join(data_folder, 'trips', '*.ekb')
    fnames = glob.glob(file_pattern)
    stmts = []
    pool = Pool(4)
    stmts_ll = pool.map(_process_trips_fname, fnames)
    pool.close()
    pool.join()
    stmts = []
    for stmts_l in stmts_ll:
        stmts += stmts_l
    return stmts


def _process_sparser_pmid(pmid):
    from indra.literature import s3_client
    try:
        logger.info('Processing %s' % pmid)
        js = s3_client.get_reader_json_str('sparser', pmid)
        jd = json.loads(js)
        sp = sparser.process_json_dict(jd)
        if sp:
            for stmt in sp.statements:
                for ev in stmt.evidence:
                    ev.pmid = pmid
            return sp.statements
        else:
            return []
    except Exception as e:
        return []


def process_sparser(data_folder):
    # Step 1: re-run reading
    from indra.tools.reading.submit_reading_pipeline import \
        submit_reading, submit_combine, wait_for_complete
    basen = 'sparser_bioexp_201905'
    job_list = submit_reading(basen, pmid_file, ['sparser'],
                              pmids_per_job=1000, force_read=True,
                              project_name='cwc')
    reading_res = wait_for_complete('run_reach_queue', job_list)
    # Step 2: re-process reading results
    from multiprocessing import Pool
    pmids = [l.strip() for l in open(pmid_file).readlines()]
    pool = Pool(4)
    stmts_ll = pool.map(_process_sparser_pmid, pmids)
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

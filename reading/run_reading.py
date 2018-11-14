import sys
import os
import time
from indra.sources.trips.drum_reader import DrumReader

if __name__ == '__main__':
    job_id = int(sys.argv[1])
    port = 6200 + job_id
    print('Starting on port %s' % port)
    pmids = [l.strip() for l in open('pmids_to_read_with_trips.txt', 'r').readlines()]
    my_pmids = pmids[job_id::50]
    for idx, pmid in enumerate(my_pmids):
        print('This is PMID #%d/%d' % (idx, len(my_pmids)))
        abs_fname = 'abstracts/%s.txt' % pmid
        if not os.path.exists(abs_fname):
            print('No abstract for %s' % pmid)
            continue
        ekb_fname = 'ekbs/%s.ekb' % pmid 
        if not os.path.exists(ekb_fname):
            print('Reading %s' % pmid)
            dr = DrumReader(run_drum=False, host='localhost', port=port)
            time.sleep(3)
            with open(abs_fname, 'r') as fh:
                abs_text = fh.read()
            if not abs_text:
                print('Empty text in abstract for %s' % pmid)
                continue
            dr.read_text(abs_text)
            try:
                dr.start()
            except SystemExit:
                if not dr.extractions:
                    print('No extractions for %s' % pmid)
                    continue
                ekb = dr.extractions[0]
                with open(ekb_fname, 'w') as fh:
                    fh.write(ekb)
            print('Finished reading %s' % pmid)
        else:
            print('Skipping existing EKB for %s' % pmid)

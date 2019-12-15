"""A script to copy the "legacy" curations from the Google Doc spreadsheet
(REACH curations, 2019-11-20) to the INDRA DB curations table."""
import sys
import csv
from indra.tools import assemble_corpus as ac
from indra_db.client.principal.curation import submit_curation
from bioexp.util import pklload

if __name__ == '__main__':
    # Get a dict of all curations by UUID
    curation_file = sys.argv[1]
    stmt_file = sys.argv[2]

    stmts = ac.load_statements(stmt_file)
    stmt_dict = {s.uuid: s for s in stmts}

    # Some default args for the curations
    source = "bioexp_paper_tsv"
    ip = '134.174.140.78'

    curations = {}
    with open(curation_file, 'r') as fh:
        reader = csv.reader(fh, delimiter='\t')
        # Skip the header row
        next(reader)
        for row in reader:
            uuid = row[1]
            ev_text = row[11]

            # Remap curator info
            curator = row[14]
            if curator == '':
                continue
            assert curator in ('BMG', 'JAB')
            curator = ('bachmanjohn@gmail.com' if curator == 'JAB'
                                               else 'ben.gyori@gmail.com')

            # Remap score field to one of the pre-set curation tags
            correct = row[15]
            if curator == '':
                continue
            elif int(correct) == 1:
                tag = 'correct'
            elif int(correct) == 0:
                tag = 'other'
            else:
                print("Invalid value for correct")
                print(row)
            comment = None if row[16] == '' else row[16]


            # Get stmt information
            stmt = stmt_dict[uuid]
            pa_hash = stmt.get_hash()

            #if pa_hash == 23109912213960991:
            #    import ipdb; ipdb.set_trace()

            stmt_ev = None
            for ev in stmt.evidence:
                if ev.text == ev_text and ev.source_api == 'reach':
                    stmt_ev = ev
            if stmt_ev is None:
                print("Could not find evidence for statement %s, text %s")
                continue
            # Something weird here--different texts are producing same
            # source hash TODO TODO TODO
            else:
                source_hash = stmt_ev.get_source_hash()

            print(pa_hash, tag, curator, ip, comment, source_hash,
                  'bioexp_paper_tsv')
            submit_curation(pa_hash, tag, curator, ip, comment, source_hash,
                            source='bioexp_paper_tsv')


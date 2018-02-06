import re
import json
from fuzzywuzzy import fuzz
from util import pklload


def fix_text(txt):
    """Eliminate some symbols to have cleaner text to read."""
    txt = re.sub('[ ]?\( xref \)', '', txt)
    # This is to make [ xref ] become [] to match the two readers
    txt = re.sub('\[ xref \]', '[]', txt)
    txt = re.sub('[\(]?XREF_BIBR[\)]?[,]?', '', txt)
    txt = re.sub('[\(]?XREF_FIG[\)]?[,]?', '', txt)
    txt = re.sub('[\(]?XREF_SUPPLEMENT[\)]?[,]?', '', txt)
    txt = txt.strip()
    return txt


if __name__ == '__main__':
    stmts = pklload('preassembled')

    text_to_read = {}
    # Iterate over all the statements
    for stmt in stmts:
        for ev in stmt.evidence:
            # Take only evidences from REACH and Sparser
            if ev.source_api in ('reach', 'sparser'):
                txt_fixed = fix_text(ev.text)
                # Create entry for PMID if it doesn't exist yet
                if ev.pmid not in text_to_read:
                    text_to_read[ev.pmid] = [txt_fixed]
                # Add evidence text to PMID entry if it's not in there yet
                elif txt_fixed not in text_to_read[ev.pmid]:
                    # Here we check whether the sentence is already in the
                    # set of sentences for the given PMID but perhaps
                    # in a trivially modified form (e.g. greek letters spelled
                    # out, hyphens missing, citation brackets added, etc.).
                    # This is done by fuzzy string matching and using a cutoff
                    # to eliminate duplicates.
                    any_similar = False
                    for txt_ref in text_to_read[ev.pmid]:
                        fuzz_ratio = fuzz.ratio(txt_ref, txt_fixed)
                        if fuzz_ratio > 90:
                            print('"%s"\n too similar to "%s"' %
                                  (txt_ref, txt_fixed))
                            any_similar = True
                            break
                    if not any_similar:
                        text_to_read[ev.pmid].append(txt_fixed)
    # Dump the sentences in a JSON file
    with open('trips_text.json', 'wt') as fh:
        json.dump(text_to_read, fh, indent=1)

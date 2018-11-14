import os
import sys
from indra.literature import s3_client


if __name__ == '__main__':
    startid = int(sys.argv[1])
    with open('pmids.txt', 'r') as fh:
        pmids = [l.strip() for l in fh.readlines()]
    for pmid in pmids[startid:startid+5000]:
        if not os.path.exists('abstracts/%s.txt' % pmid):
            print('Getting abstract %s' % pmid)
            abst = s3_client.get_gz_object('papers/PMID%s/abstract' % pmid)
            if abst is None:
                print('Failed to get abstract for %s' % pmid)
                continue
            with open('abstracts/%s.txt' % pmid, 'w') as fh:
                fh.write(abst)
        else:
            print('Skipping %s' % pmid)

"""This script counts the different content types that were read for the
Benchmark Corpus within a non-public S3-based repository."""
import json
import tqdm
import boto3
from collections import Counter, defaultdict


if __name__ == '__main__':
    s3 = boto3.client('s3')
    content_types = defaultdict(list)
    with open('pmids.txt', 'r') as fh:
        pmids = [l.strip() for l in fh.readlines()]
    irrelevants = {'reach', 'pmcid', 'doi', 'sparser', 'reach_no_bioentities'}
    for pmid in tqdm.tqdm(pmids):
        res = s3.list_objects(Bucket='bigmech', Prefix='papers/PMID%s/' % pmid)
        if 'Contents' not in res:
            content_types[pmid] = []
        else:
            contents = res['Contents']
            for content in contents:
                if any(irrelevant in content['Key']
                       for irrelevant in irrelevants):
                    continue
                elif 'abstract' in content['Key']:
                    content_types[pmid].append('abstract')
                elif 'fulltext/pmc_oa_xml' in content['Key']:
                    content_types[pmid].append('pmc_oa')
                elif 'fulltext/pmc_oa_txt' in content['Key']:
                    content_types[pmid].append('pmc_oa')
                elif 'fulltext/txt' in content['Key']:
                    content_types[pmid].append('pmc_oa')
                elif 'fulltext/nxml' in content['Key']:
                    content_types[pmid].append('pmc_oa')
                elif 'fulltext/pmc_auth_xml' in content['Key']:
                    content_types[pmid].append('pmc_manuscript')
                elif 'fulltext/elsevier_xml' in content['Key']:
                    content_types[pmid].append('elsevier')
                else:
                    print('Unhandled content key: %s' % content['Key'])

    with open('content_types_from_s3.json', 'w') as fh:
        json.dump(content_types, fh, indent=1)

    content_types_used = []
    for k, v in content_types.items():
        if not v:
            content_types_used.append('missing')
            continue
        for ct in ['pmc_oa', 'pmc_manuscript', 'elsevier', 'abstract']:
            if ct in v:
                content_types_used.append(ct)
                break

    cnt = Counter(content_types_used)
    total = sum(cnt.values())
    for k, v in cnt.most_common():
        print(k, v, '%.1f%%' % (100 * v / total))

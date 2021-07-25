import tqdm
import boto3
from collections import defaultdict


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
                elif 'fulltext/pmc_auth_xml' in content['Key']:
                    content_types[pmid].append('pmc_manuscript')
                elif 'fulltext/elsevier_xml' in content['Key']:
                    content_types[pmid].append('elsevier')
                else:
                    print('Unhandled content key: %s' % content['Key'])

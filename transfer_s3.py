import sys
import shutil
from os.path import dirname, join
import boto3

data_dir = join(dirname(__file__), 'data')

action = sys.argv[1]
if action not in ('get', 'put'):
    print("Usage: %s [get|put] <filename>" % sys.argv[0])
    sys.exit(1)

bucket_name = 'bioexp-paper'

s3_client = boto3.client('s3')

if action == 'get':
    s3_key = sys.argv[2]
    filename = join(data_dir, s3_key)
    # Get data from S3
    obj_response = s3_client.get_object(Bucket=bucket_name, Key=s3_key)
    data = obj_response['Body']
    # Write data to a file
    with open(filename, 'wb') as f:
        print("Downloading to %s..." % filename)
        shutil.copyfileobj(data, f)
else:
    filename = sys.argv[2]
    s3_key = filename
    with open(filename, 'rb') as f:
        s3_client.put_object(Bucket=bucket_name, Key=s3_key, Body=f)


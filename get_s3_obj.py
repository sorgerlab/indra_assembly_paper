from os.path import dirname, join
import sys
import boto3
import shutil

data_dir = join(dirname(__file__), 'data')

s3_key = sys.argv[1]
bucket_name = 'bioexp-paper'

s3_client = boto3.client('s3')
# Get data from S3
obj_response = s3_client.get_object(Bucket=bucket_name, Key=s3_key)
data = obj_response['Body']
# Write data to a file
filename = join(data_dir, s3_key)
with open(filename, 'wb') as f:
    print("Downloading to %s..." % filename)
    shutil.copyfileobj(data, f)

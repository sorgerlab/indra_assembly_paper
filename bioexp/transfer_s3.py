import sys
import shutil
from os.path import dirname, join, basename
import boto3


bucket_name = 'bioexp-paper'
s3_client = boto3.client('s3')


def download_from_s3(key, target_dir):
    # Get data from S3
    print("Getting %s from S3..." % key)
    obj_response = s3_client.get_object(Bucket=bucket_name, Key=key)
    data = obj_response['Body']
    filename = join(target_dir, key)
    # Write data to a file
    with open(filename, 'wb') as f:
        print("Downloading to %s..." % filename)
        shutil.copyfileobj(data, f)


def upload_to_s3(filename):
    key = basename(filename)
    print("Uploading %s to S3..." % key)
    with open(filename, 'rb') as f:
        s3_client.put_object(Bucket=bucket_name, Key=key, Body=f)


if __name__ == '__main__':
    action = sys.argv[1]
    if action not in ('get', 'put'):
        usage =  "Usage:\n"
        usage += "  %s get <key_name> <target_dir>\n" % sys.argv[0]
        usage += "  %s put <filename>" % sys.argv[0]
        print(usage)
        sys.exit(1)

    if action == 'get':
        s3_key = basename(sys.argv[2])
        target_dir = sys.argv[3]
        download_from_s3(s3_key, target_dir)
    else:
        filename = sys.argv[2]
        upload_to_s3(filename)

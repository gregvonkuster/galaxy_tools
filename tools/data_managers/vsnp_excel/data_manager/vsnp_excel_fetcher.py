import argparse
import json
import os
try:
    # For Python 3.0 and later
    from urllib.request import Request, urlopen
except ImportError:
    # Fall back to Python 2 imports
    from urllib2 import Request, urlopen


def url_download(url, target_directory):
    file_path = os.path.abspath(os.path.join(target_directory, os.path.basename(url)))
    src = None
    dst = None
    try:
        req = Request(url)
        src = urlopen(req)
        with open(file_path, 'wb') as dst:
            while True:
                chunk = src.read(2**10)
                if chunk:
                    dst.write(chunk)
                else:
                    break
    finally:
        if src:
            src.close()
    return file_path


def download(dbkey, display_name, url, out_file):

    with open(out_file) as fh:
        params = json.loads(fh.read())

    work_directory = params['output_data'][0]['extra_files_path']
    target_directory = os.path.abspath(os.path.join(work_directory, "vsnp", dbkey))
    os.mkdir(target_directory)
    file_path = url_download(url, target_directory)
    entry_name = "%s" % os.path.basename(file_path)

    data_manager_json = {"data_tables": {}}
    data_manager_entry = {}
    data_manager_entry['value'] = dbkey
    data_manager_entry['name'] = entry_name
    data_manager_entry['path'] = file_path
    data_manager_entry['description'] = "Excel file for %s" % display_name
    data_manager_json["data_tables"]["vsnp_excel"] = data_manager_entry

    with open(out_file, 'w') as fh:
        fh.write(json.dumps(data_manager_json, sort_keys=True))


parser = argparse.ArgumentParser()

parser.add_argument('--dbkey', dest='dbkey', help='Genome reference dbkey')
parser.add_argument('--display_name', dest='display_name', help='Reference display name')
parser.add_argument('--url', dest='url', help='URL to download Excel file')
parser.add_argument('--out_file', dest='out_file', help='JSON output file')

args = parser.parse_args()

download(args.dbkey, args.display_name, args.url, args.out_file)

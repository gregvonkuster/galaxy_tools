#!/usr/bin/env python
import argparse
import os
import socket
import sys
from json import dumps, loads

from six.moves.urllib.parse import urlencode, urlparse
from six.moves.urllib.request import urlopen

from galaxy.datatypes import sniff
from galaxy.datatypes.registry import Registry
from galaxy.jobs import TOOL_PROVIDED_JOB_METADATA_FILE
from galaxy.util import get_charset_from_http_headers

GALAXY_PARAM_PREFIX = "GALAXY"
GALAXY_ROOT_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
GALAXY_DATATYPES_CONF_FILE = os.path.join(GALAXY_ROOT_DIR, "datatypes_conf.xml")

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--datatypes_config", dest="datatypes_config", help="Galaxy data types conf file for mapping file types")
parser.add_argument("-l", "--latitude", dest="latitude", type=float, help="Latitude of station")
parser.add_argument("-m", "--longitude", dest="longitude", type=float, help="Latitude of station")
parser.add_argument("-o", "--output", dest="output", help="Output dataset which also acts as a temporary JSON parameter file")
parser.add_argument("-r", "--galaxy_root", dest="galaxy_root", help="Galaxy root dir")

args = parser.parse_args()

def stop_err(msg):
    sys.stderr.write(msg)
    sys.exit()


def load_input_parameters(filename, erase_file=True):
    datasource_params = {}
    try:
        json_params = loads(open(filename, "r").read())
        datasource_params = json_params.get("param_dict")
    except Exception:
        json_params = None
        for line in open(filename, "r"):
            try:
                line = line.strip()
                fields = line.split("\t")
                datasource_params[fields[0]] = fields[1]
            except Exception:
                continue
    if erase_file:
        # Ppen file for writing, then close to remove params from file.
        open(filename, "w").close()
    return json_params, datasource_params


job_params, params = load_input_parameters(args.output)
# Specially named file for output junk to pass onto set metadata.
json_file = open(job_params["job_config"]["TOOL_PROVIDED_JOB_METADATA_FILE"], "w")
datatypes_registry = Registry()
datatypes_registry.load_datatypes(root_dir=job_params["job_config"]["GALAXY_ROOT_DIR"], config=job_params["job_config"]["GALAXY_DATATYPES_CONF_FILE"])
# Using exactly URL indicates that only one dataset is being downloaded.
URL = "%s?lon=%5f&lat=%5f" % (params.get("URL", None), args.longitude, args.latitude)
URL_method = params.get("URL_method", None)
# The Python support for fetching resources from the web is layered.
# urllib uses the httplib library, which in turn uses the socket library.
# As of Python 2.3 you can specify how long a socket should wait for a
# response before timing out. By default the socket module has no timeout
# and can hang. Currently, the socket timeout is not exposed at the httplib
# or urllib2 levels. However, you can set the default timeout (in seconds)
# globally for all sockets by doing the following.
socket.setdefaulttimeout(600)

for data_dict in job_params["output_data"]:
    output_filename = data_dict.get("file_name", args.output)
    cur_URL = params.get("%s|%s|URL" % (GALAXY_PARAM_PREFIX, data_dict["out_data_name"]), URL)
    if not cur_URL or urlparse(cur_URL).scheme not in ("http", "https", "ftp"):
        open(output_filename, "w").write("")
        stop_err("The remote data source application has not sent back a URL parameter in the request.")
    # The following calls to urlopen() will use the above default timeout.
    try:
        if not URL_method or URL_method == "get":
            page = urlopen(cur_URL)
        elif URL_method == "post":
            page = urlopen(cur_URL, urlencode(params))
    except Exception as e:
        stop_err("The remote data source application may be off line, please try again later. Error: %s" % str(e))
    # Do sniff stream for multi_byte.
    try:
        output_filename, is_multi_byte = sniff.stream_to_open_named_file(page, os.open(output_filename, os.O_WRONLY | os.O_CREAT), output_filename, source_encoding=get_charset_from_http_headers(page.headers))
    except Exception as e:
        stop_err("Unable to fetch %s:\n%s" % (cur_URL, e))
    try:
        ext = sniff.handle_uploaded_dataset_file(args.output, datatypes_registry, ext=data_dict["ext"], is_multi_byte=is_multi_byte)
    except Exception as e:
        stop_err(str(e))
    info = dict(type="dataset",
                dataset_id=data_dict["dataset_id"],
                ext=ext)
    json_file.write("%s\n" % dumps(info))

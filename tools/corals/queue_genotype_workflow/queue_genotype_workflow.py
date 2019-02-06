#!/usr/bin/env python
import argparse
import json
import qgw_util

from bioblend import galaxy

parser = argparse.ArgumentParser()
parser.add_argument('--affy_metadata', dest='affy_metadata', help='Input Affymetrix 96 well plate metadata file')
parser.add_argument('--annot', dest='annot', help='Probeset annotation file')
parser.add_argument('--api_key', dest='api_key', help='Current user API key')
parser.add_argument('--calls', dest='calls', help='Apt-probeset genotype calls file')
parser.add_argument('--confidences', dest='confidences', help='Apt-probeset genotype confidences file')
parser.add_argument('--config_file', dest='config_file', help='qgw_config.ini')
parser.add_argument('--dbkey', dest='dbkey', help='Reference genome dbkey')
parser.add_argument('--reference_genome', dest='reference_genome', help='Reference genome')
parser.add_argument('--history_id', dest='history_id', help='Encoded id of current history')
parser.add_argument('--output', dest='output', help='Output dataset')
parser.add_argument('--report', dest='report', help='Apt-probeset genotype report file')
parser.add_argument('--snp-posteriors', dest='snp-posteriors', help='Apt-probeset genotype snp-posteriors file')
parser.add_argument('--summary', dest='summary', help='Apt-probeset genotype summary file')
parser.add_argument('--user_email', dest='user_email', help='Galaxy user email')
parser.add_argument('--user_id', dest='user_id', help='Encoded Galaxy user id')
args = parser.parse_args()

dh = open("/tmp/work/qgw.log", "w")
dh.write("history_id: %s\n" % str(args.history_id))
dh.write("user_id: %s\n" % str(args.user_id))
dh.write("user_email: %s\n" % str(args.user_email))
dh.write("args.api_key: %s\n" % str(args.api_key))
dh.write("args.reference_genome: %s\n" % str(args.reference_genome))
dh.write("args.dbkey: %s\n" % str(args.dbkey))
user_api_key = open(args.api_key, 'r').read()
dh.write("user_api_key: %s\n" % str(user_api_key))
galaxy_base_url = qgw_util.get_value_from_config(dh, args.config_file, 'GALAXY_BASE_URL')
dh.write("galaxy_base_url: %s\n" % str(galaxy_base_url))
gi = galaxy.GalaxyInstance(url=galaxy_base_url, key=user_api_key)
dh.write("gi: %s\n" % str(gi))
#gi_user_id = qgw_util.get_current_user(dh, gi)

workflow_name = qgw_util.get_value_from_config(dh, args.config_file, 'WORKFLOW_NAME')
dh.write("\nworkflow_name: %s\n" % str(workflow_name))

outh = open(args.output, 'w')
# Share the current history with the genotype user.
#history_dict = qgw_util.share_history(dh, gi, args.history_id, gi_user_id)

# Check the status of the last history to make
# sure it is in a finished state so that we can
# start the workflow with this history.
#histories = qgw_util.get_histories(dh, gi)
#for history in histories:
#    dh.write("history_dict:\n%s\n" % str(history))
#    history_id = history['id']
#    dh.write("history_id: %s" % str(history_id))
#    status = qgw_util.get_history_status(dh, gi, history_id)
#    dh.write("status: %s\n" % str(status))

# TODO: queue and wait until all histories are in a finished state.

# Share the history with the user.
#history_dict = qgw_util.share_history(dh, gi, history_id, args.user_id)
#outh.write("\nShared history\n%s\n with user id %s and user_email %s:\n" % (str(history_dict), args.user_id, args.user_email))

# Get the workflow.
workflow_id, workflow_dict = qgw_util.get_workflow(gi, workflow_name, dh)
dh.write("\nworkflow_id: %s\n" % str(workflow_id))
#dh.write("\nworkflow_dict: %s\n" % json.dumps(workflow_dict, sort_keys=True, indent=4))

outh.close()

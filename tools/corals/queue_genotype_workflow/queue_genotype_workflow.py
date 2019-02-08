#!/usr/bin/env python
import argparse
import json
import string
import time

from bioblend import galaxy
from six.moves import configparser

FINISHED_STATES = ['ok', 'error']

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
parser.add_argument('--history_name', dest='history_name', help='Name of current history')
parser.add_argument('--output', dest='output', help='Output dataset')
parser.add_argument('--report', dest='report', help='Apt-probeset genotype report file')
parser.add_argument('--snp-posteriors', dest='snp-posteriors', help='Apt-probeset genotype snp-posteriors file')
parser.add_argument('--summary', dest='summary', help='Apt-probeset genotype summary file')
parser.add_argument('--user_email', dest='user_email', help='Galaxy user email')
parser.add_argument('--user_id', dest='user_id', help='Encoded Galaxy user id')
args = parser.parse_args()


def add_library_dataset_to_history(gi, history_id, dataset_id, history_datasets, dh):
    # Add a data library dataset to a history.
    dh.write('\nImporting dataset id %s into history %s.\n' % (dataset_id, history_id))
    new_hda_dict = gi.histories.upload_dataset_from_library(history_id, dataset_id)
    dh.write('Response from importing dataset id %s into history %s:\n' % (dataset_id, history_id))
    dh.write('%s\n\n' % str(new_hda_dict))
    new_hda_name = new_hda_dict['name']
    history_datasets[new_hda_name] = new_hda_dict
    return history_datasets


def get_all_genotyped_samples_dataset_id(gi, data_lib_id, dh):
    """
    Use the Galaxy API to get the all samples dataset.
    We're assuming it is in the root folder.
    """
    dh.write('Searching for all_genotyped_samples dataset.\n')
    lib_item_dicts = gi.libraries.show_library(data_lib_id, contents=True)
    for lib_item_dict in lib_item_dicts:
        if lib_item_dict['type'] == 'file':
            dataset_name = lib_item_dict['name'].lstrip('/').lower()
            if dataset_name.startswith("all_genotyped_samples"):
                dh.write('Found all samples dataset.\n')
                return lib_item_dict['id']
    return None


def get_config_settings(dh, config_file, section='defaults'):
    dh.write("config_file: %s\n" % str(config_file))
    dh.write("section: %s\n" % str(section))
    d = {}
    config_parser = configparser.ConfigParser()
    config_parser.read(config_file)
    dh.write("config_parser.sections(): %s\n" % str(config_parser.sections()))
    for key, value in config_parser.items(section):
        dh.write("key: %s\n" % str(key))
        dh.write("value: %s\n" % str(value))
        if section == 'defaults':
            d[string.upper(key)] = value
        else:
            d[key] = value
    dh.write("d: %s\n" % str(d))
    return d


def get_data_library(gi, name, dh):
    # Use the Galaxy API to get the data library named the value name.
    dh.write('Searching for data library named %s.\n' % name)
    # The following is not correctly filtering out deleted libraries.
    data_lib_dicts = gi.libraries.get_libraries(library_id=None, name=name, deleted=False)
    for data_lib_dict in data_lib_dicts:
        if data_lib_dict['name'] == name and data_lib_dict['deleted'] not in [True, 'true', 'True']:
            dh.write('Found data library named %s.\n' % name)
            return data_lib_dict['id']
    return None


def get_history_datasets(gi, history_id, dh):
    dh.write('Getting datasets for history id %s\n' % history_id)
    history_datasets = {}
    history_dict = gi.histories.show_history(history_id, contents=True, deleted='false', details=None)
    #dh.write("history_dict:\n%s\n" % str(history_dict))
    for contents_dict in history_dict:
        dh.write("contents_dict:\n%s\n" % str(contents_dict))
        if contents_dict['history_content_type'] == 'dataset':
            dataset_name = contents_dict['name']
            # Don't include the "Queue genotype workflow" dataset.
            if dataset_name.startswith("Queue genotype workflow"):
                continue
            dh.write("dataset_name:\n%s\n" % str(dataset_name))
            history_datasets[dataset_name] = contents_dict
    return history_datasets


def get_value_from_config(dh, config_file, value):
    dh.write("config_file: %s\n" % str(config_file))
    dh.write("value: %s\n" % str(value))
    defaults = get_config_settings(dh, config_file)
    config_value = defaults.get(value, None)
    dh.write("config_value: %s\n" % str(config_value))
    return config_value


def get_workflow(gi, name, dh, galaxy_base_url=None, api_key=None):
    dh.write('Searching for workflow named %s\n' % name)
    workflow_info_dicts = gi.workflows.get_workflows(name=name, published=True)
    dh.write('workflow_info_dicts: %s\n' % workflow_info_dicts)
    if len(workflow_info_dicts) == 0:
        return None, None
    wf_info_dict = workflow_info_dicts[0]
    workflow_id = wf_info_dict['id']
    # Get the complete workflow.
    workflow_dict = gi.workflows.show_workflow(workflow_id)
    dh.write('Found workflow named %s.\n' % name)
    return workflow_id, workflow_dict


def get_workflow_input_datasets(gi, history_datasets, workflow_dict, dh):
    # Map the history datasets to the input datasets for the workflow.
    workflow_inputs = {}
    dh.write('\nMapping datasets from history to input datasets in workflow %s.\n' % workflow_name)
    steps_dict = workflow_dict.get('steps', None)
    dh.write("steps_dict: %s\n" % json.dumps(steps_dict, indent=4))
    if steps_dict is not None:
        for step_index, step_dict in steps_dict.items():
            # Dicts that define dataset inputs for a workflow
            # look like this.
            # "0": {
            #      "tool_id": null,
            #      "tool_version": null,
            #      "id": 0,
            #      "input_steps": {},
            #      "tool_inputs": {},
            #      "type": "data_input",
            #      "annotation": null
            # },
            tool_id = step_dict.get('tool_id', None)
            dh.write("tool_id: %s\n" % str(tool_id))
            tool_type = step_dict.get('type', None)
            dh.write("tool_type: %s\n" % str(tool_type))
            # This requires the workflow input dataset annotation to be a
            # string # (e.g., report) that enables it to be appropriatey
            # matched to a dataset (e.g., axiongt1_report.txt).
            # 1. affy_metadata.tabular - must have the word "metadata" in
            #                            the annotation.
            # 2. sample_attributes.txt - must have the word "attributes" in
            #                            the annotation.
            # 3. probeset_annotation.csv - must have the word "annotation" in
            #                              the annotation.
            # 4. <summary file>.txt - must have the the word "summary" in the
            #                         annotation.
            # 5. <snp-posteriors file>.txt - must have the the word
            #                                "snp-posteriors" in the annotation.
            # 6. <report file>.txt - must have the the word "report" in the
            #                        annotation.
            # 7. <confidences file>.txt - must have the the word "confidences"
            #                             in the annotation.
            # 8. <calls file>.txt - must have the the word "calls" in the
            #                       annotation.
            # TODO: figure out how to name the "all samples" dataset so
            # we can properly map it here.
            annotation = step_dict.get('annotation', None)
            dh.write("annotation: %s\n" % str(annotation))
            if tool_id is None and tool_type == 'data_input' and annotation is not None:
                annotation_check = annotation.lower()
                # inputs is a list and workflow input datasets
                # have no inputs.
                for input_hda_name, input_hda_dict in history_datasets.items():
                    input_hda_name_check = input_hda_name.lower()
                    dh.write("input_hda_name_check: %s\n" % str(input_hda_name_check))
                    if input_hda_name_check.find(annotation_check) >= 0:
                        workflow_inputs[step_index] = {'src': 'hda', 'id': input_hda_dict['id']}
                        dh.write('Mapped dataset %s from history to workflow input dataset with annotation %s.\n' % (input_hda_name, annotation))
                        break
    return workflow_inputs


def start_workflow(gi, workflow_id, workflow_name, inputs, params, history_id, dh):
    dh.write('\nExecuting workflow %s.\n' % workflow_name)
    dh.write('inputs:\n%s\n\n' % str(inputs))
    dh.write('history_id:\n%s\n\n' % str(history_id))
    workflow_invocation_dict = gi.workflows.invoke_workflow(workflow_id, inputs=inputs, params=params, history_id=history_id)
    dh.write('Response from executing workflow %s:\n' % workflow_name)
    dh.write('%s\n' % str(workflow_invocation_dict))


def update_workflow_params(workflow_dict, dbkey, dh):
    parameter_updates = None
    name = workflow_dict['name']
    dh.write('Checking for tool parameter updates for workflow %s using dbkey %s.\n' % (name, dbkey))
    step_dicts = workflow_dict.get('steps', None)
    for step_id, step_dict in step_dicts.items():
        tool_id = step_dict['tool_id']
        if tool_id is None:
            continue
        # Handle reference_source entries
        if tool_id.find('affy2vcf') > 0:
            dh.write('\nChecking tool id: %s\n' % tool_id)
            tool_inputs_dict = step_dict['tool_inputs']
            # The queue_genotype_workflow tool provides a selection of only
            # a locally cached reference genome (not a history item), so dbkey
            # will always refer to a locally cached genome.
            # The affy2vcf tool allows the user to select either a locally
            # cached reference genome or a history item, but the workflow is
            # defined to use a locally cached reference genome by default.
            reference_genome_source_cond_dict = tool_inputs_dict['reference_genome_source_cond']
            # The value of reference_genome_source_cond_dict['reference_genome_source']
            # will always be 'cached'.
            workflow_db_key = reference_genome_source_cond_dict['locally_cached_item']
            if dbkey !=  workflow_db_key:
                reference_genome_source_cond_dict['locally_cached_item'] = dbkey
                parameter_updates = {}
                parameter_updates[step_id] = reference_genome_source_cond_dict
                dh.write('Updated step id %s with the following entry:\n%s\n' % (step_id, str(reference_genome_source_cond_dict)))
    return parameter_updates


dh = open("/tmp/work/qgw.log", "w")
dh.write("history_id: %s\n" % str(args.history_id))
dh.write("history_name: %s\n" % str(args.history_name))
dh.write("user_id: %s\n" % str(args.user_id))
dh.write("user_email: %s\n" % str(args.user_email))
dh.write("args.api_key: %s\n" % str(args.api_key))
dh.write("args.reference_genome: %s\n" % str(args.reference_genome))
dh.write("args.dbkey: %s\n" % str(args.dbkey))
user_api_key = open(args.api_key, 'r').read()
dh.write("user_api_key: %s\n" % str(user_api_key))
galaxy_base_url = get_value_from_config(dh, args.config_file, 'GALAXY_BASE_URL')
dh.write("galaxy_base_url: %s\n" % str(galaxy_base_url))
gi = galaxy.GalaxyInstance(url=galaxy_base_url, key=user_api_key)
dh.write("gi: %s\n" % str(gi))
# gi_user_id = get_current_user(dh, gi)

all_genotyped_samples_library_name = get_value_from_config(dh, args.config_file, 'ALL_GENOTYPED_SAMPLES_LIBRARY_NAME')
workflow_name = get_value_from_config(dh, args.config_file, 'WORKFLOW_NAME')
#dh.write("\nworkflow_name: %s\n" % str(workflow_name))

outh = open(args.output, 'w')

# Get the workflow.
workflow_id, workflow_dict = get_workflow(gi, workflow_name, dh)
#dh.write("\nworkflow_id: %s\n" % str(workflow_id))
# dh.write("\nworkflow_dict: %s\n" % json.dumps(workflow_dict, sort_keys=True, indent=4))

# Get the All Genotyped Samples data library.
all_genotyped_samples_library_id = get_data_library(gi, all_genotyped_samples_library_name, dh)
#dh.write("\nall_genotyped_samples_library_id: %s\n" % str(all_genotyped_samples_library_id))

# Get the public all_genotyped_samples" dataset id.
all_genotyped_samples_dataset_id = get_all_genotyped_samples_dataset_id(gi, all_genotyped_samples_library_id, dh)
#dh.write("\nall_genotyped_samples_dataset_id: %s\n" % str(all_genotyped_samples_dataset_id))

# Get the current history datasets.  At this point, the history must contain
# only the datasets to be used as input to the workflow.
history_datasets = get_history_datasets(gi, args.history_id, dh)
#dh.write("\nhistory_datasets: %s\n" % str(history_datasets))
# dh.write("\nSleeping for 10 seconds...\n")
# time.sleep(10)

# Import the public all_genotyped_samples dataset from the data library to the history.
# TODO uncomment the following when we're ready.  it defintiely works!
# history_datasets = add_library_dataset_to_history(gi, args.history_id, all_genotyped_samples_dataset_id, history_datasets, dh)
# dh.write("\nhistory_datasets: %s\n" % str(history_datasets))

# Map the history datasets to the input datasets for the workflow.
workflow_input_datasets = get_workflow_input_datasets(gi, history_datasets, workflow_dict, dh)
#dh.write("\nworkflow_input_datasets: %s\n" % str(workflow_input_datasets))
# dh.write("\nSleeping for 10 seconds...\n")
# time.sleep(10)

# Get teh workflow params that could be updated.
params = update_workflow_params(workflow_dict, args.dbkey, dh) 
dh.write("\nparams:\n%s\n" % str(params))

# Start the workflow.
start_workflow(gi, workflow_id, workflow_name, workflow_input_datasets, params, args.history_id, dh)
dh.write("\nSleeping for 10 seconds...\n")
time.sleep(10)

# TODO: the following won't work because the history state is "running" as
# long as this tool is not finished.  So we should probably add a tool to
# the end of the workflow that removes a file that is created by this tool.
# Then the following loop can wait unti that file no longer exists and then
# break.

# Poll the history status and wait until it is in a finished state.
# while True:
#    status = get_history_status(dh, gi, args.history_id)
#    dh.write("status: %s\n" % str(status))
#    if status in FINISHED_STATES:
#        break
#    time.sleep(60)

outh.close()

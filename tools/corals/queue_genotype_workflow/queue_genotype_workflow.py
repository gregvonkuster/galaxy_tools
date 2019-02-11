#!/usr/bin/env python
import argparse
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
parser.add_argument('--output', dest='output', help='Output dataset')
parser.add_argument('--report', dest='report', help='Apt-probeset genotype report file')
parser.add_argument('--snp-posteriors', dest='snp-posteriors', help='Apt-probeset genotype snp-posteriors file')
parser.add_argument('--summary', dest='summary', help='Apt-probeset genotype summary file')
parser.add_argument('--user_id', dest='user_id', help='Encoded Galaxy user id')
args = parser.parse_args()


def add_library_dataset_to_history(gi, history_id, dataset_id, history_datasets, outputfh):
    # Add a data library dataset to a history.
    outputfh.write('\nImporting dataset into current history.\n')
    new_hda_dict = gi.histories.upload_dataset_from_library(history_id, dataset_id)
    new_hda_name = new_hda_dict['name']
    history_datasets[new_hda_name] = new_hda_dict
    return history_datasets


def get_all_genotyped_samples_dataset_id(gi, data_lib_id, outputfh):
    """
    Use the Galaxy API to get the all samples dataset.
    We're assuming it is in the root folder.
    """
    outputfh.write('\nSearching for dataset named all_genotyped_samples.\n')
    lib_item_dicts = gi.libraries.show_library(data_lib_id, contents=True)
    for lib_item_dict in lib_item_dicts:
        if lib_item_dict['type'] == 'file':
            dataset_name = lib_item_dict['name'].lstrip('/').lower()
            if dataset_name.startswith("all_genotyped_samples"):
                outputfh.write('Found dataset named all_genotyped_samples.\n')
                return lib_item_dict['id']
    return None


def get_config_settings(config_file, section='defaults'):
    d = {}
    config_parser = configparser.ConfigParser()
    config_parser.read(config_file)
    for key, value in config_parser.items(section):
        if section == 'defaults':
            d[string.upper(key)] = value
        else:
            d[key] = value
    return d


def get_data_library(gi, name, outputfh):
    # Use the Galaxy API to get the data library named the value name.
    outputfh.write('\nSearching for data library named %s.\n' % name)
    # The following is not correctly filtering out deleted libraries.
    data_lib_dicts = gi.libraries.get_libraries(library_id=None, name=name, deleted=False)
    for data_lib_dict in data_lib_dicts:
        if data_lib_dict['name'] == name and data_lib_dict['deleted'] not in [True, 'true', 'True']:
            outputfh.write('Found data library named %s.\n' % name)
            return data_lib_dict['id']
    return None


def get_history_status(gi, history_id):
    return gi.histories.get_status(history_id)


def get_history_datasets(gi, history_id):
    history_datasets = {}
    history_dict = gi.histories.show_history(history_id, contents=True, deleted='false', details=None)
    for contents_dict in history_dict:
        if contents_dict['history_content_type'] == 'dataset':
            dataset_name = contents_dict['name']
            # Don't include the "Queue genotype workflow" dataset.
            if dataset_name.startswith("Queue genotype workflow"):
                continue
            history_datasets[dataset_name] = contents_dict
    return history_datasets


def get_value_from_config(config_file, value):
    defaults = get_config_settings(config_file)
    config_value = defaults.get(value, None)
    return config_value


def get_workflow(gi, name, outputfh, galaxy_base_url=None, api_key=None):
    outputfh.write('\nSearching for workflow named %s\n' % name)
    workflow_info_dicts = gi.workflows.get_workflows(name=name, published=True)
    if len(workflow_info_dicts) == 0:
        return None, None
    wf_info_dict = workflow_info_dicts[0]
    workflow_id = wf_info_dict['id']
    # Get the complete workflow.
    workflow_dict = gi.workflows.show_workflow(workflow_id)
    outputfh.write('Found workflow named %s.\n' % name)
    return workflow_id, workflow_dict


def get_workflow_input_datasets(gi, history_datasets, workflow_dict, outputfh):
    # Map the history datasets to the input datasets for the workflow.
    workflow_inputs = {}
    outputfh.write('\nMapping datasets from history to workflow %s.\n' % workflow_name)
    steps_dict = workflow_dict.get('steps', None)
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
            tool_type = step_dict.get('type', None)
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
            if tool_id is None and tool_type == 'data_input' and annotation is not None:
                annotation_check = annotation.lower()
                # inputs is a list and workflow input datasets
                # have no inputs.
                for input_hda_name, input_hda_dict in history_datasets.items():
                    input_hda_name_check = input_hda_name.lower()
                    if input_hda_name_check.find(annotation_check) >= 0:
                        workflow_inputs[step_index] = {'src': 'hda', 'id': input_hda_dict['id']}
                        outputfh.write(' - Mapped dataset %s from history to workflow input dataset with annotation %s.\n' % (input_hda_name, annotation))
                        break
    return workflow_inputs


def start_workflow(gi, workflow_id, workflow_name, inputs, params, history_id, outputfh):
    outputfh.write('\nExecuting workflow %s.\n' % workflow_name)
    workflow_invocation_dict = gi.workflows.invoke_workflow(workflow_id, inputs=inputs, params=params, history_id=history_id)
    outputfh.write('Response from executing workflow %s:\n' % workflow_name)
    outputfh.write('%s\n' % str(workflow_invocation_dict))


def update_workflow_params(workflow_dict, dbkey, outputfh):
    parameter_updates = None
    name = workflow_dict['name']
    outputfh.write('\nChecking for tool parameter updates for workflow %s using dbkey %s.\n' % (name, dbkey))
    step_dicts = workflow_dict.get('steps', None)
    for step_id, step_dict in step_dicts.items():
        tool_id = step_dict['tool_id']
        if tool_id is None:
            continue
        # Handle reference_source entries
        if tool_id.find('affy2vcf') > 0:
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
            if dbkey != workflow_db_key:
                reference_genome_source_cond_dict['locally_cached_item'] = dbkey
                parameter_updates = {}
                parameter_updates[step_id] = reference_genome_source_cond_dict
                outputfh.write('Updated step id %s with the following entry:\n%s\n' % (step_id, str(reference_genome_source_cond_dict)))
    return parameter_updates


outputfh = open(args.output, "w")
user_api_key = open(args.api_key, 'r').read()
galaxy_base_url = get_value_from_config(args.config_file, 'GALAXY_BASE_URL')
gi = galaxy.GalaxyInstance(url=galaxy_base_url, key=user_api_key)
all_genotyped_samples_library_name = get_value_from_config(args.config_file, 'ALL_GENOTYPED_SAMPLES_LIBRARY_NAME')
workflow_name = get_value_from_config(args.config_file, 'WORKFLOW_NAME')

# Get the workflow.
workflow_id, workflow_dict = get_workflow(gi, workflow_name, outputfh)
# Get the All Genotyped Samples data library.
all_genotyped_samples_library_id = get_data_library(gi, all_genotyped_samples_library_name, outputfh)
# Get the public all_genotyped_samples" dataset id.
all_genotyped_samples_dataset_id = get_all_genotyped_samples_dataset_id(gi, all_genotyped_samples_library_id, outputfh)
# Get the current history datasets.  At this point, the history must contain
# only the datasets to be used as input to the workflow.
history_datasets = get_history_datasets(gi, args.history_id)
# Import the public all_genotyped_samples dataset from the data library to the history.
history_datasets = add_library_dataset_to_history(gi, args.history_id, all_genotyped_samples_dataset_id, history_datasets, outputfh)
# Map the history datasets to the input datasets for the workflow.
workflow_input_datasets = get_workflow_input_datasets(gi, history_datasets, workflow_dict, outputfh)
# Get the workflow params that could be updated.
params = update_workflow_params(workflow_dict, args.dbkey, outputfh)
# Start the workflow.
start_workflow(gi, workflow_id, workflow_name, workflow_input_datasets, params, args.history_id, outputfh)
time.sleep(10)
# Poll the history datasets, checking the statuses, and wait until
# the workflow is finished.  The workflow itself simply schedules
# all of the jobs, so it cannot be checked for a state.
while True:
    history_status_dict = get_history_status(gi, args.history_id)
    state_details_dict = history_status_dict['state_details']
    # The queue_genotype_workflow tool will continue to be in a
    # "running" state while inside this for loop, so  we know that
    # the workflow has completed if no datasets are in the "new" or
    # "queued" state and there is only 1 dataset in the "running"
    # state.  We cannot filter on datasets in the "paused" state
    # because any datasets downstream from one in an "error" state
    # will automatically be given a "paused" state.
    if state_details_dict['queued'] == 0 and state_details_dict['new'] == 0 and state_details_dict['running'] <= 1:
        break
    time.sleep(5)

outputfh.close()

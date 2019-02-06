import fileinput
import json
import numpy
import os
import shlex
import string
import subprocess
import sys
import tempfile

from six.moves import configparser
from six.moves.urllib.error import HTTPError, URLError
from six.moves.urllib.request import Request, urlopen
from six import string_types
from bioblend import galaxy


CONFIG_FILE = 'qgw_config.ini'
# Allows characters that are escaped to be un-escaped.
MAPPED_CHARS = {'>': '__gt__',
                '<': '__lt__',
                "'": '__sq__',
                '"': '__dq__',
                '[': '__ob__',
                ']': '__cb__',
                '{': '__oc__',
                '}': '__cc__',
                '@': '__at__',
                '\n': '__cn__',
                '\r': '__cr__',
                '\t': '__tc__',
                '#': '__pd__'}
# Maximum value of a signed 32 bit integer (2**31 - 1).
MAX_GENOME_SIZE = 2147483647


def add_library_dataset_to_history(gi, dbkey, history_id, history_name, dataset_id, history_input_datasets, dh):
    """
    Add a data library dataset to a history.
    """
    dh.write('\nImporting dataset id %s for dbkey %s into history %s.\n' % (dataset_id, dbkey, history_name))
    new_hda_dict = gi.histories.upload_dataset_from_library(history_id, dataset_id)
    dh.write('Response from importing dataset id %s for dbkey %s into history %s:\n' % (dataset_id, dbkey, history_name))
    dh.write('%s\n\n' % str(new_hda_dict))
    new_hda_name = new_hda_dict['name']
    history_input_datasets[new_hda_name] = new_hda_dict
    return history_input_datasets


def check_response(pegr_url, payload, response):
    try:
        s = json.dumps(payload)
        response_code = response.get('response_code', None)
        if response_code not in ['200']:
            err_msg = 'Error sending statistics to PEGR!\n\nPEGR URL:\n%s\n\n' % str(pegr_url)
            err_msg += 'Payload:\n%s\n\nResponse:\n%s\n' % (s, str(response))
            if response_code in ['500']:
                # The payload may not have included all items
                # required by PEGR, so write the error but
                # don't exit.
                sys.stderr.write(err_msg)
            else:
                # PEGR is likely unavailable, so exit.
                stop_err(err_msg)
    except Exception as e:
        err_msg = 'Error handling response from PEGR!\n\nException:\n%s\n\n' % str(e)
        err_msg += 'PEGR URL:\n%s\n\nPayload:\n%s\n\nResponse:\n%s\n' % (pegr_url, s, str(response))
        sys.stderr.write(err_msg)


def check_samtools():
    samtools_exec = which('samtools')
    if not samtools_exec:
        stop_err('Attempting to use functionality requiring samtools, but it cannot be located on Galaxy\'s PATH.')


def create_history(dh, gi, user_email, outh):
    # Create a new history to contain the analysis
    history_name = get_new_history_name(dh, gi, user_email)
    new_history_dict = gi.histories.create_history(name=history_name)
    new_history_id = new_history_dict['id']
    outh.write('\nCreated a new history named %s to contain the analysis.\n' % history_name)
    return history_name, new_history_id


def format_tool_parameters(parameters):
    s = parameters.lstrip('__SeP__')
    items = s.split('__SeP__')
    params = {}
    param_index = 0
    for i in range(len(items) / 2):
        params[restore_text(items[param_index])] = restore_text(items[param_index + 1])
        param_index += 2
    return params


def get_base_json_dict(config_file, dbkey, history_id, history_name, stats_tool_id, stderr, tool_id, tool_parameters, user_email, workflow_step_id):
    d = {}
    d['genome'] = dbkey
    d['historyId'] = history_id
    d['parameters'] = format_tool_parameters(tool_parameters)
    # d['run'] = get_run_from_history_name(history_name)
    # d['sample'] = get_sample_from_history_name(history_name)
    d['statsToolId'] = stats_tool_id
    d['toolCategory'] = get_tool_category(config_file, tool_id)
    d['toolStderr'] = stderr
    d['toolId'] = tool_id
    d['userEmail'] = user_email
    d['workflowId'] = get_workflow_id(config_file, history_name)
    d['workflowStepId'] = workflow_step_id
    return d


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


def get_current_user(dh, gi):
    current_user = gi.users.get_current_user()
    dh.write("Current user email: %s\n" % str(current_user['email']))
    return current_user['id']

def get_datasets(config_file, ids, datatypes):
    # http://localhost:8763/datasets/eca0af6fb47bf90c/display/?preview=True
    defaults = get_config_settings(config_file=config_file)
    d = {}
    for i, t in zip(listify(ids), listify(datatypes)):
        d['id'] = i
        d['type'] = t
        d['uri'] = '%s/datasets/%s/display?preview=True' % (defaults['GALAXY_BASE_URL'], i)
    return d


def get_galaxy_instance(api_key, url):
    return galaxy.GalaxyInstance(url=url, key=api_key)


def get_galaxy_url(config_file):
    defaults = get_config_settings(config_file=config_file)
    return make_url(defaults['GALAXY_API_KEY'], defaults['GALAXY_BASE_URL'])


def get_histories(dh, gi):
    histories = gi.histories.get_histories()
    dh.write("len(histories): %d\n" % len(histories))
    return histories


def get_history_index(dh, gi):
    # All history names  for the genotype user
    # are auto-generated using the format:
    # <integer index>_<user email>
    histories = get_histories(dh, gi)
    # Set a default value.
    index = 0
    for history in histories:
        name = history.get("name", None)
        dh.write("name: %s\n" % str(name))
        if name is not None:
            items = name.split("_")
            dh.write("items: %s\n" % str(items))
            if len(items) == 2:
                try:
                    index = int(items[0])
                    dh.write("index: %s\n" % str(index))
                    break
                except Exception as e:
                    # pass
                    dh.write("Exception: %s\n" % e)
    # Increment the index.
    index += 1
    return str(index).zfill(6)


def get_history_status(dh, gi, history_id):
    status = gi.histories.get_status(history_id)
    dh.write("status: %s\n" % str(status))
    return status


def get_new_history_name(dh, gi, user_email):
    dh.write("user_email: %s\n" % str(user_email))
    history_index = get_history_index(dh, gi)
    dh.write("history_index: %s\n" % str(history_index))
    history_name = "%s_%s" % (history_index, user_email)
    dh.write("history_name: %s\n" % str(history_name))
    return history_name


def get_number_of_lines(file_path):
    i = 0
    with open(file_path) as fh:
        for i, l in enumerate(fh):
            pass
    fh.close()
    if i == 0:
        return i
    return i + 1


def get_tmp_filename(dir=None, suffix=None):
    fd, name = tempfile.mkstemp(suffix=suffix, dir=dir)
    os.close(fd)
    return name


def get_tool_category(config_file, tool_id):
    lc_tool_id = tool_id.lower()
    category_map = get_config_settings(config_file=config_file)
    return category_map.get(lc_tool_id, 'Unknown')


def get_value_from_config(dh, config_file, value):
    dh.write("config_file: %s\n" % str(config_file))
    dh.write("value: %s\n" % str(value))
    defaults = get_config_settings(dh, config_file)
    config_value = defaults.get(value, None)
    dh.write("config_value: %s\n" % str(config_value))
    return config_value


def get_workflow(gi, name, dh, galaxy_base_url=None, api_key=None, for_inputs=False):
    dh.write('Searching for workflow named %s\n' % name)
    workflow_info_dicts = gi.workflows.get_workflows(name=name, published=True)
    dh.write('workflow_info_dicts: %s\n' % workflow_info_dicts)
    if len(workflow_info_dicts) == 0:
        return None, None
    wf_info_dict = workflow_info_dicts[0]
    workflow_id = wf_info_dict['id']
    # Get the complete workflow.
    if for_inputs:
        # Bioblend does not provides this end
        # point so we need to get it from Galaxy.
        base_url = '%s/api/workflows/%s/download' % (galaxy_base_url, workflow_id)
        url = api_util.make_url(api_key, base_url)
        workflow_dict = api_util.get(url)
    else:
        workflow_dict = gi.workflows.show_workflow(workflow_id)
    dh.write('Found workflow named %s.\n' % name)
    return workflow_id, workflow_dict


def get_workflow_id(config_file, history_name):
    workflow_name = get_workflow_name_from_history_name(history_name)
    if workflow_name == 'unknown':
        return 'unknown'
    defaults = get_config_settings(config_file=config_file)
    gi = get_galaxy_instance(defaults['GALAXY_API_KEY'], defaults['GALAXY_BASE_URL'])
    workflow_info_dicts = gi.workflows.get_workflows(name=workflow_name)
    if len(workflow_info_dicts) == 0:
        return 'unknown'
    wf_info_dict = workflow_info_dicts[0]
    return wf_info_dict['id']


def get_workflow_input_datasets(gi, history_name, history_input_datasets, workflow_name, dbkey, galaxy_base_url, api_key, dh):
    # Map the history datasets to the input datasets for the workflow.
    workflow_id, workflow_dict = get_workflow(gi, workflow_name, dh, galaxy_base_url=galaxy_base_url, api_key=api_key, for_inputs=True)
    workflow_inputs = {}
    dh.write('\nMapping datasets from history %s to input datasets in workflow %s.\n' % (history_name, workflow_name))
    steps_dict = workflow_dict.get('steps', None)
    if steps_dict is not None:
        for step_index, step_dict in steps_dict.items():
            inputs = step_dict.get('inputs', None)
            if inputs is not None and len(inputs) == 0:
                # inputs is a list and workflow input datasets
                # have no inputs.
                label = step_dict.get('label', None)
                if label is not None:
                    for input_hda_name, input_hda_dict in history_input_datasets.items():
                        # This requires the workflow input dataset label to be a string
                        # (e.g., R1) that is contained in the name of the input dataset
                        # (e.g., 60642_R1.fq).  The blacklist filter dataset must have
                        # the exact label "blacklist" (without the quotes).
                        if input_hda_name.find(label) >= 0:
                            workflow_inputs[step_index] = {'src': 'hda', 'id': input_hda_dict['id']}
                            dh.write('Mapped dataset %s from history to workflow input dataset with label %s.\n' % (input_hda_name, label))
                            break
    return workflow_inputs


def get_workflow_name_from_history_name(history_name, exit_on_error=False):
    # Example: paired_001-199-10749.001
    items = history_name.split('-')
    try:
        workflow_name = items[0]
    except Exception as e:
        if exit_on_error:
            stop_err('History name is likely invalid, it does not contain a workflow name: %s' % str(e))
        return 'unknown'
    return workflow_name


def listify(item, do_strip=False):
    """
    Make a single item a single item list, or return a list if passed a
    list.  Passing a None returns an empty list.
    """
    if not item:
        return []
    elif isinstance(item, list):
        return item
    elif isinstance(item, string_types) and item.count(','):
        if do_strip:
            return [token.strip() for token in item.split(',')]
        else:
            return item.split(',')
    else:
        return [item]


def make_url(api_key, url, args=None):
    """
    Adds the API Key to the URL if it's not already there.
    """
    if args is None:
        args = []
    argsep = '&'
    if '?' not in url:
        argsep = '?'
    if '?apiKey=' not in url and '&apiKey=' not in url:
        args.insert(0, ('apiKey', api_key))
    return url + argsep + '&'.join(['='.join(t) for t in args])


def post(api_key, url, data):
    url = make_url(api_key, url)
    response = Request(url, headers={'Content-Type': 'application/json'}, data=json.dumps(data))
    return json.loads(urlopen(response).read())


def restore_text(text, character_map=MAPPED_CHARS):
    """Restores sanitized text"""
    if not text:
        return text
    for key, value in character_map.items():
        text = text.replace(value, key)
    return text


def share_history(dh, gi, history_id, user_id):
    # Share a history with a user.
    history_dict = gi.histories.share_history(history_id, user_id)
    dh.write('\nShared history with id %s with the user with id %s.\n' % (history_id, user_id))
    return history_dict


def start_workflow(gi, workflow_id, workflow_name, inputs, params, history_id, dh):
    dh.write('\nExecuting workflow %s.\n' % workflow_name)
    dh.write('inputs:\n%s\n\n' % str(inputs))
    dh.write('params:\n%s\n\n' % str(params))
    dh.write('history_id:\n%s\n\n' % str(history_id))
    workflow_invocation_dict = gi.workflows.invoke_workflow(workflow_id,
                                                            inputs=inputs,
                                                            params=params,
                                                            history_id=history_id)
    dh.write('Response from executing workflow %s:\n' % workflow_name)
    dh.write('%s\n' % str(workflow_invocation_dict))


def stop_err(msg):
    sys.stderr.write(msg)
    sys.exit()


def store_results(file_path, pegr_url, payload, response):
    with open(file_path, 'w') as fh:
        # Eliminate the API key from the PEGR url.
        items = pegr_url.split('?')
        fh.write("pegr_url:\n%s\n\n" % str(items[0]))
        fh.write("payload:\n%s\n\n" % json.dumps(payload))
        fh.write("response:\n%s\n" % str(response))
        fh.close()


def submit(config_file, data):
    """
    Sends an API POST request and acts as a generic formatter for the JSON response.
    'data' will become the JSON payload read by Galaxy.
    """
    defaults = get_config_settings(config_file=config_file)
    try:
        return post(defaults['PEGR_API_KEY'], defaults['PEGR_URL'], data)
    except HTTPError as e:
        return json.loads(e.read())
    except URLError as e:
        return dict(response_code=None, message=str(e))
    except Exception as e:
        try:
            return dict(response_code=None, message=e.read())
        except:
            return dict(response_code=None, message=str(e))


def update_dataset(gi, dbkey, history_id, history_name, history_input_datasets, dh):
    """
    Update information about a history dataset.
    """
    for hda_name, attributes in history_input_datasets.items():
        if hda_name.find('blacklist') < 0:
            dataset_id = attributes['id']
            dh.write('\nUpdating dataset id %s for dbkey %s in history %s.\n' % (dataset_id, dbkey, history_name))
            # The response here is an HTTP response
            # (e.g., 200) when it should be a JSON dict.
            response = gi.histories.update_dataset(history_id=history_id, dataset_id=dataset_id, genome_build=dbkey)
            dh.write('Response from updating dataset id %s for dbkey %s in history %s:\n' % (dataset_id, dbkey, history_name))
            dh.write('%s\n\n' % str(response))
            # Since the call to update_dataset above unfortunately
            # returns an HTTP code, we need to request the updated
            # JSON object for the hda.
            new_hda_dict = gi.histories.show_dataset(history_id=history_id, dataset_id=dataset_id)
            dh.write('JSON dict from updating dataset id %s for dbkey %s in history %s:\n' % (dataset_id, dbkey, history_name))
            dh.write('%s\n\n' % str(new_hda_dict))
            new_hda_name = new_hda_dict['name']
            history_input_datasets[new_hda_name] = new_hda_dict
    return history_input_datasets


def which(file):
    # http://stackoverflow.com/questions/5226958/which-equivalent-function-in-python
    for path in os.environ["PATH"].split(":"):
        if os.path.exists(path + "/" + file):
            return path + "/" + file
    return None

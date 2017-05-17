import logging


log = logging.getLogger(__name__)

DEFAULT_DESTINATION = 'local_multicore'


def dynamic_processors_select(app, job, resource_params):
    destination = app.job_config.get_destination(DEFAULT_DESTINATION)
    if 'processors' in resource_params:
        if resource_params['processors'] < 16:
            resource_params['processors'] = 16
        destination.params['nativeSpecification'] = '--nodes=1 --ntasks=%i' % resource_params['processors']
        log.debug("dynamic_processors_select returning destination '%s' with nativeSpecification '%s'", DEFAULT_DESTINATION, destination.params['nativeSpecification'])
    else:
        log.warning("dynamic_process_select called but 'processors' not found in resource_params, returning destination '%s' with default nativeSpecification '%s'", DEFAULT_DESTINATION, destination.params['nativeSpecification'])
    return destination

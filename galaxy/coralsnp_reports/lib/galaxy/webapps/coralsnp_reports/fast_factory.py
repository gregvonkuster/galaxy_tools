"""Module containing factory class for building uvicorn app for the Galaxy CoralSNP Reports.

Information on uvicorn, its various settings, and how to invoke it can
be found at https://www.uvicorn.org/.

The Galaxy CoralSNP Reports can be launched with uvicorn using the following invocation:

::

    uvicorn --app-dir lib --factory galaxy.webapps.coralsnp_reports.fast_factory:factory

Use the environment variable ``GALAXY_CORALSNP_REPORTS_CONFIG`` to specify a Galaxy CoralSNP Reports
configuration file.

::

    GALAXY_CORALSNP_REPORTS_CONFIG=config/coralsnp_reports.yml uvicorn --app-dir lib --factory galaxy.webapps.coralsnp_reports.fast_factory:factory

.. note::

    Information on additional ways to configure uvicorn can be found at
    https://www.uvicorn.org/.


`Gunicorn <https://docs.gunicorn.org/en/stable/index.html>`__ is a server with
more complex management options.

This factory function can be executed as a uvicorn worker managed with gunicorn
with the following command-line.

::

    gunicorn 'galaxy.webapps.coralsnp_reports.fast_factory:factory()' --env GALAXY_CORALSNP_REPORTS_CONFIG=config/coralsnp_reports.yml --pythonpath lib -w 4 -k uvicorn.workers.UvicornWorker

"""

from galaxy.main_config import (
    WebappConfigResolver,
    WebappSetupProps,
)
from galaxy.webapps.coralsnp_reports.buildapp import app_factory
from .fast_app import initialize_fast_app


def factory():
    props = WebappSetupProps(
        app_name="coralsnp_reports",
        default_section_name="coralsnp_reports",
        env_config_file="GALAXY_CORALSNP_REPORTS_CONFIG",
    )
    config_provider = WebappConfigResolver(props)
    config = config_provider.resolve_config()
    gx_webapp = app_factory(
        global_conf=config.global_conf, load_app_kwds=config.load_app_kwds, wsgi_preflight=config.wsgi_preflight
    )
    return initialize_fast_app(gx_webapp)

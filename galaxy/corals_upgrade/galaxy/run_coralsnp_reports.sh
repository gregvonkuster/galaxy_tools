#!/bin/sh


# Usage: ./run_coralsnp_reports.sh [--sync-config] <start|stop>
#
#
# Description: This script can be used to start or stop the galaxy
# coralsnp reports web application. Passing in --sync-config as the first
# argument to this will cause Galaxy's database and path parameters
# from galaxy.ini to be copied over into reports.ini.

cd "$(dirname "$0")"

GALAXY_CORALSNP_REPORTS_PID=${GALAXY_CORALSNP_REPORTS_PID:-coralsnp_reports_webapp.pid}
GALAXY_CORALSNP_REPORTS_LOG=${GALAXY_CORALSNP_REPORTS_LOG:-coralsnp_reports_webapp.log}
PID_FILE=$GALAXY_CORALSNP_REPORTS_PID
LOG_FILE=$GALAXY_CORALSNP_REPORTS_LOG

. ./scripts/common_startup_functions.sh

parse_common_args $@

run_common_start_up

setup_python

if [ -z "$GALAXY_CORALSNP_REPORTS_CONFIG" ]; then
    GALAXY_CORALSNP_REPORTS_CONFIG=$(PYTHONPATH=lib python -c "from __future__ import print_function; from galaxy.util.properties import find_config_file; print(find_config_file(['coralsnp_reports', 'coralsnp_reports_wsgi']) or '')")
    export GALAXY_REPORTS_CONFIG
fi

find_server ${GALAXY_CORALSNP_REPORTS_CONFIG:-none} coralsnp_reports

if [ "$run_server" = "gunicorn" -a -z "$GALAXY_CORALSNP_REPORTS_CONFIG" ]; then
    GALAXY_CORALSNP_REPORTS_CONFIG="config/coralsnp_reports.yml.sample"
    export GALAXY_CORALSNP_REPORTS_CONFIG
    echo 'WARNING: Using default coralsnp reports config at config/coralsnp_reports.yml.sample, copy to config/coralsnp_reports.yml or set $GALAXY_CORALSNP_REPORTS_CONFIG if this is not intentional'
fi

echo "Executing: $run_server $server_args"
eval $run_server $server_args

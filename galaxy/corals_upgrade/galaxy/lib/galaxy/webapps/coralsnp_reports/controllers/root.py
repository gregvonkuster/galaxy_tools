import logging

from galaxy.webapps.base.controller import (
    BaseUIController,
    web,
)
from galaxy.webapps.coralsnp_reports.controllers.query import ReportQueryBuilder

log = logging.getLogger(__name__)


class CoralSNPReport(BaseUIController, ReportQueryBuilder):
    @web.expose
    def index(self, trans, **kwd):
        return trans.fill_template('/webapps/coralsnp_reports/index.mako')

    @web.expose
    def home(self, trans, **kwd):
        return trans.fill_template('/webapps/coralsnp_reports/home.mako')

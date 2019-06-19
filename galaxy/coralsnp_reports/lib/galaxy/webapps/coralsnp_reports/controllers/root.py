import logging

from . import BaseUIController
from galaxy.web.base.controller import web

log = logging.getLogger(__name__)


class CoralSNPReport(BaseUIController):
    @web.expose
    def index(self, trans, **kwd):
        return trans.fill_template('/webapps/coralsnp_reports/index.mako')

    @web.expose
    def home(self, trans, **kwd):
        return trans.fill_template('/webapps/coralsnp_reports/home.mako')

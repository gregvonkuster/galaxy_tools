import logging

import sqlalchemy as sa
from markupsafe import escape

import galaxy.model
from galaxy import util
from . import BaseUIController
from galaxy.web.base.controller import web
from galaxy.webapps.reports.controllers.query import ReportQueryBuilder

log = logging.getLogger(__name__)


class Collectors(BaseUIController, ReportQueryBuilder):

    @web.expose
    def all(self, trans, **kwd):
        message = escape(util.restore_text(kwd.get('message', '')))
        q = sa.select((galaxy.model.corals.Person.table.c.id,
                       galaxy.model.corals.Person.table.c.last_name,
                       galaxy.model.corals.Person.table.c.first_name,
                       galaxy.model.corals.Person.table.c.organization,
                       galaxy.model.corals.Person.table.c.email),
                      from_obj=[galaxy.model.corals.Person.table],
                      order_by=[galaxy.model.corals.Person.table.c.last_name])
        collectors = []
        for row in q.execute():
            cols_tup = (row.id, row.last_name, row.first_name, row.organization, row.email)
            collectors.append(cols_tup)
        return trans.fill_template('/webapps/coralsnp_reports/collectors.mako', collectors=collectors, message=message)

    @web.expose
    def of_sample(self, trans, **kwd):
        message = escape(util.restore_text(kwd.get('message', '')))
        affy_id = kwd.get('affy_id')
        collector_id = kwd.get('collector_id')
        q = sa.select((galaxy.model.corals.Person.table.c.id,
                       galaxy.model.corals.Person.table.c.last_name,
                       galaxy.model.corals.Person.table.c.first_name,
                       galaxy.model.corals.Person.table.c.organization,
                       galaxy.model.corals.Person.table.c.email),
                      whereclause=galaxy.model.corals.Person.table.c.id == collector_id,
                      from_obj=[galaxy.model.corals.Person.table])
        collectors = []
        for row in q.execute():
            cols_tup = (row.last_name, row.first_name, row.organization, row.email)
            collectors.append(cols_tup)
        return trans.fill_template('/webapps/coralsnp_reports/collector_of_sample.mako',
                                   affy_id=affy_id,
                                   collectors=collectors,
                                   message=message)

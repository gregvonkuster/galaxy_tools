import calendar
import logging
import operator
from datetime import (
    date,
    datetime,
    timedelta
)

import sqlalchemy as sa
from markupsafe import escape
from sqlalchemy import false

import galaxy.model
from galaxy import util
from . import BaseUIController
from galaxy.web.base.controller import web
from galaxy.webapps.reports.controllers.jobs import sorter
from galaxy.webapps.reports.controllers.query import ReportQueryBuilder

log = logging.getLogger(__name__)


class Colonies(BaseUIController, ReportQueryBuilder):

    @web.expose
    def all(self, trans, **kwd):
        message = escape(util.restore_text(kwd.get('message', '')))
        q = sa.select((galaxy.model.corals.Colony.table.c.id,
                       galaxy.model.corals.Colony.table.c.latitude,
                       galaxy.model.corals.Colony.table.c.longitude,
                       galaxy.model.corals.Colony.table.c.depth,
                       galaxy.model.corals.Colony.table.c.reef_id),
                      from_obj=[galaxy.model.corals.Colony.table],
                      order_by=[galaxy.model.corals.Colony.table.c.id])
        colonies = []
        for row in q.execute():
            cols_tup = (row.id, row.latitude, row.longitude,
                        row.depth, row.reef_id)
            colonies.append(cols_tup)
        return trans.fill_template('/webapps/coralsnp_reports/colonies.mako', colonies=colonies, message=message)

    @web.expose
    def of_sample(self, trans, **kwd):
        message = escape(util.restore_text(kwd.get('message', '')))
        affy_id = kwd.get('affy_id')
        colony_id = kwd.get('colony_id')
        q = sa.select((galaxy.model.corals.Colony.table.c.latitude,
                       galaxy.model.corals.Colony.table.c.longitude,
                       galaxy.model.corals.Colony.table.c.depth,
                       galaxy.model.corals.Colony.table.c.reef_id),
                      whereclause=sa.and_(galaxy.model.corals.Colony.table.c.id == colony_id),
                      from_obj=[galaxy.model.corals.Colony.table],
                      order_by=[galaxy.model.corals.Colony.table.c.id])
        colonies = []
        for row in q.execute():
            cols_tup = (row.latitude, row.longitude, row.depth, row.reef_id)
            colonies.append(cols_tup)
        return trans.fill_template('/webapps/coralsnp_reports/colony_of_sample.mako',
                                   affy_id=affy_id,
                                   colonies=colonies,
                                   message=message)

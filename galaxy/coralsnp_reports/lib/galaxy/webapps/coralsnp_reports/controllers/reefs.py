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


class Reefs(BaseUIController, ReportQueryBuilder):

    @web.expose
    def all(self, trans, **kwd):
        message = escape(util.restore_text(kwd.get('message', '')))
        q = sa.select((galaxy.model.corals.Reef.table.c.id,
                       galaxy.model.corals.Reef.table.c.name,
                       galaxy.model.corals.Reef.table.c.region,
                       galaxy.model.corals.Reef.table.c.latitude,
                       galaxy.model.corals.Reef.table.c.longitude,
                       galaxy.model.corals.Reef.table.c.geographic_origin),
                      from_obj=[galaxy.model.corals.Reef.table],
                      order_by=[galaxy.model.corals.Reef.table.c.id])
        reefs = []
        for row in q.execute():
            cols_tup = (row.id, row.name, row.region, row.latitude, row.longitude, row.geographic_origin)
            reefs.append(cols_tup)
        return trans.fill_template('/webapps/coralsnp_reports/reefs.mako', reefs=reefs, message=message)

    @web.expose
    def of_colony(self, trans, **kwd):
        message = escape(util.restore_text(kwd.get('message', '')))
        reef_id = kwd.get('reef_id')
        latitude = kwd.get('latitude')
        longitude = kwd.get('longitude')
        depth = kwd.get('depth')
        q = sa.select((galaxy.model.corals.Reef.table.c.name,
                       galaxy.model.corals.Reef.table.c.region,
                       galaxy.model.corals.Reef.table.c.latitude,
                       galaxy.model.corals.Reef.table.c.longitude,
                       galaxy.model.corals.Reef.table.c.geographic_origin),
                      whereclause=sa.and_(galaxy.model.corals.Reef.table.c.id == reef_id),
                      from_obj=[galaxy.model.corals.Reef.table],
                      order_by=[galaxy.model.corals.Reef.table.c.id])
        reefs = []
        for row in q.execute():
            cols_tup = (row.name, row.region, row.latitude, row.longitude, row.geographic_origin)
            reefs.append(cols_tup)
        return trans.fill_template('/webapps/coralsnp_reports/reef_of_colony.mako',
                                   latitude=latitude,
                                   longitude=longitude,
                                   depth=depth,
                                   reefs=reefs,
                                   message=message)

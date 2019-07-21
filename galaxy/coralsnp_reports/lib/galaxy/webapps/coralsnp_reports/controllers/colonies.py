import logging

import sqlalchemy as sa
from markupsafe import escape

import galaxy.model
from galaxy import util
from . import BaseUIController
from galaxy.web.base.controller import web
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
                       galaxy.model.corals.Colony.table.c.reef_id,
                       galaxy.model.corals.Sample.table.c.public,
                       galaxy.model.corals.Sample.table.c.public_after_date),
                      from_obj=[galaxy.model.corals.Colony.table,
                                galaxy.model.corals.Sample.table],
                      whereclause=galaxy.model.corals.Colony.table.c.id == galaxy.model.corals.Sample.table.c.colony_id,
                      order_by=[galaxy.model.corals.Colony.table.c.id])
        colonies = []
        for row in q.execute():
            public_after_date = str(row.public_after_date)[:10]
            if str(row.public) == "True":
                latitude = row.latitude
                longitude = row.longitude
            else:
                latitude = "Private until %s" % public_after_date
                longitude = "Private until %s" % public_after_date
            cols_tup = (row.id, latitude, longitude, row.depth, row.reef_id)
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
                       galaxy.model.corals.Colony.table.c.reef_id,
                       galaxy.model.corals.Sample.table.c.public,
                       galaxy.model.corals.Sample.table.c.public_after_date),
                      from_obj=[galaxy.model.corals.Colony.table,
                                galaxy.model.corals.Sample.table],
                      whereclause=sa.and_(galaxy.model.corals.Colony.table.c.id == colony_id,
                                          galaxy.model.corals.Colony.table.c.id == galaxy.model.corals.Sample.table.c.colony_id),
                      order_by=[galaxy.model.corals.Colony.table.c.id])
        colonies = []
        for row in q.execute():
            public_after_date = str(row.public_after_date)[:10]
            if str(row.public) == "True":
                latitude = row.latitude
                longitude = row.longitude
            else:
                latitude = "Private until %s" % public_after_date
                longitude = "Private until %s" % public_after_date
            cols_tup = (latitude, longitude, row.depth, row.reef_id)
            colonies.append(cols_tup)
        return trans.fill_template('/webapps/coralsnp_reports/colony_of_sample.mako',
                                   affy_id=affy_id,
                                   colonies=colonies,
                                   message=message)

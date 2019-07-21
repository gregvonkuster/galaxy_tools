import logging

import sqlalchemy as sa
from markupsafe import escape

import galaxy.model
from galaxy import util
from . import BaseUIController
from galaxy.web.base.controller import web
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
                       galaxy.model.corals.Reef.table.c.geographic_origin,
                       galaxy.model.corals.Sample.table.c.colony_id,
                       galaxy.model.corals.Sample.table.c.public,
                       galaxy.model.corals.Sample.table.c.public_after_date),
                      from_obj=[galaxy.model.corals.Reef.table,
                                galaxy.model.corals.Colony.table,
                                galaxy.model.corals.Sample.table],
                      whereclause=sa.and_(galaxy.model.corals.Sample.table.c.colony_id == galaxy.model.corals.Colony.table.c.id,
                                          galaxy.model.corals.Colony.table.c.reef_id == galaxy.model.corals.Reef.table.c.id),
                      order_by=[galaxy.model.corals.Reef.table.c.id])
        reefs = []
        for row in q.execute():
            public_after_date = str(row.public_after_date)[:10]
            if str(row.public) == "True":
                latitude = row.latitude
                longitude = row.longitude
            else:
                latitude = "Private until %s" % public_after_date
                longitude = "Private until %s" % public_after_date
            cols_tup = (row.id, row.name, row.region, latitude, longitude, row.geographic_origin)
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
                       galaxy.model.corals.Reef.table.c.geographic_origin,
                       galaxy.model.corals.Sample.table.c.public,
                       galaxy.model.corals.Sample.table.c.public_after_date),
                      from_obj=[galaxy.model.corals.Reef.table,
                                galaxy.model.corals.Colony.table,
                                galaxy.model.corals.Sample.table],
                      whereclause=sa.and_(galaxy.model.corals.Reef.table.c.id == reef_id,
                                          galaxy.model.corals.Sample.table.c.colony_id == galaxy.model.corals.Colony.table.c.id,
                                          galaxy.model.corals.Colony.table.c.reef_id == galaxy.model.corals.Reef.table.c.id),
                      order_by=[galaxy.model.corals.Reef.table.c.id])
        reefs = []
        for row in q.execute():
            public_after_date = str(row.public_after_date)[:10]
            if str(row.public) == "True":
                latitude = row.latitude
                longitude = row.longitude
            else:
                latitude = "Private until %s" % public_after_date
                longitude = "Private until %s" % public_after_date
            cols_tup = (row.name, row.region, latitude, longitude, row.geographic_origin)
            reefs.append(cols_tup)
        return trans.fill_template('/webapps/coralsnp_reports/reef_of_colony.mako',
                                   latitude=latitude,
                                   longitude=longitude,
                                   depth=depth,
                                   reefs=reefs,
                                   message=message)

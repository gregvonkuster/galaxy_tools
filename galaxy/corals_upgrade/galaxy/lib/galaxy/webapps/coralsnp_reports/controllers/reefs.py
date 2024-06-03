import logging

import sqlalchemy as sa
from markupsafe import escape

from galaxy.model import corals
from galaxy import util
from galaxy.webapps.base.controller import (
    BaseUIController,
    web,
)
from galaxy.webapps.coralsnp_reports.controllers.query import ReportQueryBuilder

log = logging.getLogger(__name__)


class Reefs(BaseUIController, ReportQueryBuilder):

    @web.expose
    def all(self, trans, **kwd):
        message = escape(util.restore_text(kwd.get('message', '')))
        q = sa.select(
            corals.Reef.id,
            corals.Reef.name,
            corals.Reef.region,
            corals.Reef.latitude,
            corals.Reef.longitude,
            corals.Reef.geographic_origin,
            corals.Sample.colony_id,
            corals.Sample.public,
            corals.Sample.public_after_date
        ).select_from(corals.Reef).join(corals.Colony, corals.Reef.id == corals.Colony.reef_id).join(corals.Sample, corals.Colony.id == corals.Sample.colony_id).order_by(corals.Reef.id)
        reefs = []
        for row in trans.sa_session.execute(q):
            public_after_date = str(row.public_after_date)[:10]
            if str(row.public) == "True":
                latitude = row.latitude
                longitude = row.longitude
            else:
                latitude = "Private until %s" % public_after_date
                longitude = "Private until %s" % public_after_date
            reefs.append((row.id, row.name, row.region, latitude, longitude, row.geographic_origin))
        return trans.fill_template('/webapps/coralsnp_reports/reefs.mako', reefs=reefs, message=message)

    @web.expose
    def of_colony(self, trans, **kwd):
        message = escape(util.restore_text(kwd.get('message', '')))
        reef_id = kwd.get('reef_id')
        latitude = kwd.get('latitude')
        longitude = kwd.get('longitude')
        depth = kwd.get('depth')
        q = sa.select(
            corals.Reef.name,
            corals.Reef.region,
            corals.Reef.latitude,
            corals.Reef.longitude,
            corals.Reef.geographic_origin,
            corals.Sample.public,
            corals.Sample.public_after_date
        ).select_from(corals.Reef).where(corals.Reef.id == reef_id).join(corals.Colony, corals.Reef.id == corals.Colony.reef_id).join(corals.Sample, corals.Sample.colony_id == corals.Colony.id).order_by(corals.Reef.id)
        reefs = []
        for row in trans.sa_session.execute(q):
            public_after_date = str(row.public_after_date)[:10]
            if str(row.public) == "True":
                latitude = row.latitude
                longitude = row.longitude
            else:
                latitude = "Private until %s" % public_after_date
                longitude = "Private until %s" % public_after_date
            reefs.append((row.name, row.region, latitude, longitude, row.geographic_origin))
        return trans.fill_template('/webapps/coralsnp_reports/reef_of_colony.mako',
                                   latitude=latitude,
                                   longitude=longitude,
                                   depth=depth,
                                   reefs=reefs,
                                   message=message)

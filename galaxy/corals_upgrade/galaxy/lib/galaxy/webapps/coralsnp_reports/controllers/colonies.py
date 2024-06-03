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


class Colonies(BaseUIController, ReportQueryBuilder):

    @web.expose
    def all(self, trans, **kwd):
        message = escape(util.restore_text(kwd.get('message', '')))
        q = sa.select(
            corals.Colony.id,
            corals.Colony.latitude,
            corals.Colony.longitude,
            corals.Colony.depth,
            corals.Colony.reef_id,
            corals.Sample.public,
            corals.Sample.public_after_date
        ).select_from(corals.Colony).join(corals.Sample, corals.Colony.id == corals.Sample.colony_id).order_by(corals.Colony.id)
        colonies = []
        for row in trans.sa_session.execute(q):
            public_after_date = str(row.public_after_date)[:10]
            if str(row.public) == "True":
                latitude = row.latitude
                longitude = row.longitude
            else:
                latitude = "Private until %s" % public_after_date
                longitude = "Private until %s" % public_after_date
            colonies.append((row.id, latitude, longitude, row.depth, row.reef_id))
        return trans.fill_template('/webapps/coralsnp_reports/colonies.mako',
                                   colonies=colonies,
                                   message=message)

    @web.expose
    def of_sample(self, trans, **kwd):
        message = escape(util.restore_text(kwd.get('message', '')))
        affy_id = kwd.get('affy_id')
        colony_id = kwd.get('colony_id')
        q = sa.select(
            corals.Colony.latitude,
            corals.Colony.longitude,
            corals.Colony.depth,
            corals.Colony.reef_id,
            corals.Sample.public,
            corals.Sample.public_after_date).select_from(corals.Colony).join(corals.Sample, sa.and_(corals.Colony.id == colony_id, corals.Colony.id == corals.Sample.colony_id)).order_by(corals.Colony.id)
        colonies = []
        for row in trans.sa_session.execute(q):
            public_after_date = str(row.public_after_date)[:10]
            if str(row.public) == "True":
                latitude = row.latitude
                longitude = row.longitude
            else:
                latitude = "Private until %s" % public_after_date
                longitude = "Private until %s" % public_after_date
            colonies.append((latitude, longitude, row.depth, row.reef_id))
        return trans.fill_template('/webapps/coralsnp_reports/colony_of_sample.mako',
                                   affy_id=affy_id,
                                   colonies=colonies,
                                   message=message)

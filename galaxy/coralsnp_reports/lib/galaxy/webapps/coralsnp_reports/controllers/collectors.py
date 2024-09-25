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


class Collectors(BaseUIController, ReportQueryBuilder):

    @web.expose
    def all(self, trans, **kwd):
        message = escape(util.restore_text(kwd.get('message', '')))
        collectors = []
        for row in trans.sa_session.query(corals.Person):
            collectors.append((row.id, row.last_name, row.first_name, row.organization, row.email))
        return trans.fill_template('/webapps/coralsnp_reports/collectors.mako',
                                   collectors=collectors,
                                   message=message)

    @web.expose
    def of_sample(self, trans, **kwd):
        message = escape(util.restore_text(kwd.get('message', '')))
        affy_id = kwd.get('affy_id')
        collector_id = kwd.get('collector_id')
        q = (
            sa.select(
                corals.Person.id,
                corals.Person.last_name,
                corals.Person.first_name,
                corals.Person.organization,
                corals.Person.email
            )
            .select_from(corals.Person.table)
            .where(corals.Person.table.c.id == collector_id)
        )
        collectors = []
        for row in trans.sa_session.execute(q):
            collectors.append((row.last_name, row.first_name, row.organization, row.email))
        return trans.fill_template('/webapps/coralsnp_reports/collector_of_sample.mako',
                                   affy_id=affy_id,
                                   collectors=collectors,
                                   message=message)

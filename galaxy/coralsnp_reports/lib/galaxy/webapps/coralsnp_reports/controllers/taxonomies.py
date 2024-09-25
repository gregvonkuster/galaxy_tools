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


class Taxonomies(BaseUIController, ReportQueryBuilder):

    @web.expose
    def all(self, trans, **kwd):
        message = escape(util.restore_text(kwd.get('message', '')))
        taxonomies = []
        for row in trans.sa_session.query(corals.Taxonomy):
            taxonomies.append((row.id, row.species_name, row.genus_name))
        return trans.fill_template('/webapps/coralsnp_reports/taxonomies.mako', taxonomies=taxonomies, message=message)

    @web.expose
    def of_sample(self, trans, **kwd):
        message = escape(util.restore_text(kwd.get('message', '')))
        affy_id = kwd.get('affy_id')
        taxonomy_id = kwd.get('taxonomy_id')
        q = (
            sa.select(
                corals.Taxonomy.species_name,
                corals.Taxonomy.genus_name
            )
            .select_from(corals.Taxonomy.table)
            .where(corals.Taxonomy.table.c.id == taxonomy_id)
            .order_by(corals.Taxonomy.table.c.id)
        )
        taxonomies = []
        for row in trans.sa_session.execute(q):
            taxonomies.append((row.species_name, row.genus_name))
        return trans.fill_template('/webapps/coralsnp_reports/taxonomy_of_sample.mako',
                                   affy_id=affy_id,
                                   taxonomies=taxonomies,
                                   message=message)

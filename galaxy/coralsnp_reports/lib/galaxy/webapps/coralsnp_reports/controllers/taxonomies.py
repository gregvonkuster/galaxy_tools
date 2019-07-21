import logging

import sqlalchemy as sa
from markupsafe import escape

import galaxy.model
from galaxy import util
from . import BaseUIController
from galaxy.web.base.controller import web
from galaxy.webapps.reports.controllers.query import ReportQueryBuilder

log = logging.getLogger(__name__)


class Taxonomies(BaseUIController, ReportQueryBuilder):

    @web.expose
    def all(self, trans, **kwd):
        message = escape(util.restore_text(kwd.get('message', '')))
        q = sa.select((galaxy.model.corals.Taxonomy.table.c.id,
                       galaxy.model.corals.Taxonomy.table.c.species_name,
                       galaxy.model.corals.Taxonomy.table.c.genus_name),
                      from_obj=[galaxy.model.corals.Taxonomy.table],
                      order_by=[galaxy.model.corals.Taxonomy.table.c.id])
        taxonomies = []
        for row in q.execute():
            cols_tup = (row.id, row.species_name, row.genus_name)
            taxonomies.append(cols_tup)
        return trans.fill_template('/webapps/coralsnp_reports/taxonomies.mako', taxonomies=taxonomies, message=message)

    @web.expose
    def of_sample(self, trans, **kwd):
        message = escape(util.restore_text(kwd.get('message', '')))
        affy_id = kwd.get('affy_id')
        taxonomy_id = kwd.get('taxonomy_id')
        q = sa.select((galaxy.model.corals.Taxonomy.table.c.species_name,
                       galaxy.model.corals.Taxonomy.table.c.genus_name),
                      whereclause=sa.and_(galaxy.model.corals.Taxonomy.table.c.id == taxonomy_id),
                      from_obj=[galaxy.model.corals.Taxonomy.table],
                      order_by=[galaxy.model.corals.Taxonomy.table.c.id])
        taxonomies = []
        for row in q.execute():
            cols_tup = (row.species_name, row.genus_name)
            taxonomies.append(cols_tup)
        return trans.fill_template('/webapps/coralsnp_reports/taxonomy_of_sample.mako',
                                   affy_id=affy_id,
                                   taxonomies=taxonomies,
                                   message=message)

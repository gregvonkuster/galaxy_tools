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


class Genotypes(BaseUIController, ReportQueryBuilder):

    @web.expose
    def all(self, trans, **kwd):
        message = escape(util.restore_text(kwd.get('message', '')))
        num_genotypes = trans.sa_session.query(galaxy.model.corals.Genotype).count()
        return trans.fill_template('/webapps/coralsnp_reports/genotypes.mako', num_genotypes=num_genotypes, message=message)

    @web.expose
    def for_sample(self, trans, **kwd):
        message = escape(util.restore_text(kwd.get('message', '')))
        # If specified_date is not received, we'll default to the current month
        affy_id = kwd.get('affy_id')
        genotype_id = kwd.get('genotype_id')
        sample_id = kwd.get('sample_id')
        user_specimen_id = kwd.get('user_specimen_id')
        q = sa.select((galaxy.model.corals.Genotype.table.c.coral_mlg_clonal_id,
                       galaxy.model.corals.Genotype.table.c.coral_mlg_rep_sample_id,
                       galaxy.model.corals.Genotype.table.c.genetic_coral_species_call,
                       galaxy.model.corals.Genotype.table.c.bcoral_genet_id),
                      whereclause=sa.and_(galaxy.model.corals.Genotype.table.c.id == genotype_id),
                      from_obj=[galaxy.model.corals.Genotype.table],
                      order_by=[galaxy.model.corals.Genotype.table.c.id])
        genotypes = []
        for row in q.execute():
            cols_tup = (row.coral_mlg_clonal_id, row.coral_mlg_rep_sample_id,
                        row.genetic_coral_species_call, row.bcoral_genet_id)
            genotypes.append((cols_tup))
        return trans.fill_template('/webapps/coralsnp_reports/genotypes_for_sample.mako',
                                   affy_id=affy_id,
                                   genotypes=genotypes,
                                   sample_id=sample_id,
                                   user_specimen_id=user_specimen_id,
                                   message=message)

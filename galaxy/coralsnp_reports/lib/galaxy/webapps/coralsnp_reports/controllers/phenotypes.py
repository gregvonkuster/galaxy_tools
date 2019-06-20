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


class Phenotypes(BaseUIController, ReportQueryBuilder):

    @web.expose
    def all(self, trans, **kwd):
        message = escape(util.restore_text(kwd.get('message', '')))
        q = sa.select((galaxy.model.corals.Phenotype.table.c.id,
                       galaxy.model.corals.Phenotype.table.c.disease_resist,
                       galaxy.model.corals.Phenotype.table.c.bleach_resist,
                       galaxy.model.corals.Phenotype.table.c.mortality,
                       galaxy.model.corals.Phenotype.table.c.tle,
                       galaxy.model.corals.Phenotype.table.c.spawning,
                       galaxy.model.corals.Phenotype.table.c.sperm_motility,
                       galaxy.model.corals.Phenotype.table.c.healing_time),
                      from_obj=[galaxy.model.corals.Phenotype.table],
                      order_by=[galaxy.model.corals.Phenotype.table.c.id])
        phenotypes = []
        for row in q.execute():
            cols_tup = (row.id, row.disease_resist, row.bleach_resist,
                        row.mortality, row.tle, row.spawning, row.sperm_motility,
                        row.healing_time)
            phenotypes.append(cols_tup)
        return trans.fill_template('/webapps/coralsnp_reports/phenotypes.mako', phenotypes=phenotypes, message=message)

    @web.expose
    def for_sample(self, trans, **kwd):
        message = escape(util.restore_text(kwd.get('message', '')))
        affy_id = kwd.get('affy_id')
        phenotype_id = kwd.get('phenotype_id')
        q = sa.select((galaxy.model.corals.Phenotype.table.c.disease_resist,
                       galaxy.model.corals.Phenotype.table.c.bleach_resist,
                       galaxy.model.corals.Phenotype.table.c.mortality,
                       galaxy.model.corals.Phenotype.table.c.tle,
                       galaxy.model.corals.Phenotype.table.c.spawning,
                       galaxy.model.corals.Phenotype.table.c.sperm_motility,
                       galaxy.model.corals.Phenotype.table.c.healing_time),
                      whereclause=sa.and_(galaxy.model.corals.Phenotype.table.c.id == phenotype_id),
                      from_obj=[galaxy.model.corals.Phenotype.table],
                      order_by=[galaxy.model.corals.Phenotype.table.c.id])
        phenotypes = []
        for row in q.execute():
            cols_tup = (row.disease_resist, row.bleach_resist, row.mortality, row.tle,
                        row.spawning, row.sperm_motility, row.healing_time)
            phenotypes.append(cols_tup)
        return trans.fill_template('/webapps/coralsnp_reports/phenotype_of_sample.mako',
                                   affy_id=affy_id,
                                   phenotypes=phenotypes,
                                   message=message)

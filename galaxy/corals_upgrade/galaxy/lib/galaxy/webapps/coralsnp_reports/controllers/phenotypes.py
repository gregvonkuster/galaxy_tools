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


class Phenotypes(BaseUIController, ReportQueryBuilder):

    @web.expose
    def all(self, trans, **kwd):
        message = escape(util.restore_text(kwd.get('message', '')))
        phenotypes = []
        for row in trans.sa_session.query(corals.Phenotype):
            phenotypes.append((row.id, row.disease_resist, row.bleach_resist, row.mortality, row.tle, row.spawning, row.sperm_motility, row.healing_time))
        return trans.fill_template('/webapps/coralsnp_reports/phenotypes.mako', phenotypes=phenotypes, message=message)

    @web.expose
    def of_sample(self, trans, **kwd):
        message = escape(util.restore_text(kwd.get('message', '')))
        affy_id = kwd.get('affy_id')
        phenotype_id = kwd.get('phenotype_id')
        q = sa.select(
            (corals.Phenotype.disease_resist,
             corals.Phenotype.bleach_resist,
             corals.Phenotype.mortality,
             corals.Phenotype.tle,
             corals.Phenotype.spawning,
             corals.Phenotype.sperm_motility,
             corals.Phenotype.healing_time),
            from_obj=[corals.Phenotype.table],
            whereclause=(corals.Phenotype.table.c.id == phenotype_id),
            order_by=[corals.Phenotype.table.c.id])
        phenotypes = []
        for row in trans.sa_session.execute(q):
            phenotypes.append((row.disease_resist, row.bleach_resist, row.mortality, row.tle, row.spawning, row.sperm_motility, row.healing_time))
        return trans.fill_template('/webapps/coralsnp_reports/phenotype_of_sample.mako',
                                   affy_id=affy_id,
                                   phenotypes=phenotypes,
                                   message=message)

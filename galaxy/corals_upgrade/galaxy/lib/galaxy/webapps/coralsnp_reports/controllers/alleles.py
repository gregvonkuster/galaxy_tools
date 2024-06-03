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


class Alleles(BaseUIController, ReportQueryBuilder):

    @web.expose
    def of_sample(self, trans, **kwd):
        message = escape(util.restore_text(kwd.get('message', '')))
        affy_id = kwd.get('affy_id')
        allele_id = kwd.get('allele_id')
        q = sa.select(
            (corals.Allele.id,
             corals.Allele.allele),
            from_obj=[corals.Allele.table],
            whereclause=corals.Allele.table.c.id == allele_id,
        )
        alleles = []
        for row in trans.sa_session.execute(q):
            alleles.append((row.id, row.allele))
        return trans.fill_template('/webapps/coralsnp_reports/allele_of_sample.mako',
                                   affy_id=affy_id,
                                   alleles=alleles,
                                   message=message)

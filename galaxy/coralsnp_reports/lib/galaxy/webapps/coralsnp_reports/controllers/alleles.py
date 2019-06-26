import logging

import sqlalchemy as sa
from markupsafe import escape

import galaxy.model
from galaxy import util
from . import BaseUIController
from galaxy.web.base.controller import web
from galaxy.webapps.reports.controllers.query import ReportQueryBuilder

log = logging.getLogger(__name__)


class Alleles(BaseUIController, ReportQueryBuilder):

    @web.expose
    def of_sample(self, trans, **kwd):
        message = escape(util.restore_text(kwd.get('message', '')))
        affy_id = kwd.get('affy_id')
        allele_id = kwd.get('allele_id')
        q = sa.select((galaxy.model.corals.Allele.table.c.id,
                       galaxy.model.corals.Allele.table.c.allele),
                      whereclause=galaxy.model.corals.Allele.table.c.id == allele_id,
                      from_obj=[galaxy.model.corals.Allele.table])
        alleles = []
        for row in q.execute():
            cols_tup = (row.id, row.allele)
            alleles.append(cols_tup)
        return trans.fill_template('/webapps/coralsnp_reports/allele_of_sample.mako',
                                   affy_id=affy_id,
                                   alleles=alleles,
                                   message=message)

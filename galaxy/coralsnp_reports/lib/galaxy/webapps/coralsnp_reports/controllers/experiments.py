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


class Experiments(BaseUIController, ReportQueryBuilder):

    @web.expose
    def all(self, trans, **kwd):
        message = escape(util.restore_text(kwd.get('message', '')))
        q = sa.select((galaxy.model.corals.Experiment.table.c.id,
                       galaxy.model.corals.Experiment.table.c.seq_facility,
                       galaxy.model.corals.Experiment.table.c.array_version,
                       galaxy.model.corals.Experiment.table.c.result_folder_name,
                       galaxy.model.corals.Experiment.table.c.plate_barcode),
                      from_obj=[galaxy.model.corals.Experiment.table],
                      order_by=[galaxy.model.corals.Experiment.table.c.id])
        experiments = []
        for row in q.execute():
            cols_tup = (row.id, row.seq_facility, row.array_version,
                        row.result_folder_name, row.plate_barcode)
            experiments.append(cols_tup)
        return trans.fill_template('/webapps/coralsnp_reports/experiments.mako', experiments=experiments, message=message)

    @web.expose
    def of_sample(self, trans, **kwd):
        message = escape(util.restore_text(kwd.get('message', '')))
        affy_id = kwd.get('affy_id')
        experiment_id = kwd.get('experiment_id')
        q = sa.select((galaxy.model.corals.Experiment.table.c.seq_facility,
                       galaxy.model.corals.Experiment.table.c.array_version,
                       galaxy.model.corals.Experiment.table.c.result_folder_name,
                       galaxy.model.corals.Experiment.table.c.plate_barcode),
                      whereclause=sa.and_(galaxy.model.corals.Experiment.table.c.id == experiment_id),
                      from_obj=[galaxy.model.corals.Experiment.table],
                      order_by=[galaxy.model.corals.Experiment.table.c.id])
        experiments = []
        for row in q.execute():
            cols_tup = (row.seq_facility, row.array_version,
                        row.result_folder_name, row.plate_barcode)
            experiments.append(cols_tup)
        return trans.fill_template('/webapps/coralsnp_reports/experiment_of_sample.mako',
                                   affy_id=affy_id,
                                   experiments=experiments,
                                   message=message)

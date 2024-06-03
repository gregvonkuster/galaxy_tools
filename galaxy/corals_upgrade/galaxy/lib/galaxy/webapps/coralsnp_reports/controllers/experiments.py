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


class Experiments(BaseUIController, ReportQueryBuilder):

    @web.expose
    def all(self, trans, **kwd):
        message = escape(util.restore_text(kwd.get('message', '')))
        experiments = []
        for row in trans.sa_session.query(corals.Experiment):
            experiments.append((row.id, row.seq_facility, row.array_version, row.result_folder_name, row.plate_barcode))
        return trans.fill_template('/webapps/coralsnp_reports/experiments.mako', experiments=experiments, message=message)

    @web.expose
    def of_sample(self, trans, **kwd):
        message = escape(util.restore_text(kwd.get('message', '')))
        affy_id = kwd.get('affy_id')
        experiment_id = kwd.get('experiment_id')
        q = sa.select(
            (corals.Experiment.seq_facility,
             corals.Experiment.array_version,
             corals.Experiment.result_folder_name,
             corals.Experiment.plate_barcode),
            from_obj=[corals.Experiment.table],
            whereclause=corals.Experiment.table.c.id == experiment_id,
            order_by=[corals.Experiment.table.c.id]
        )
        experiments = []
        for row in trans.sa_session.execute(q):
            experiments.append((row.seq_facility, row.array_version, row.result_folder_name, row.plate_barcode))
        return trans.fill_template('/webapps/coralsnp_reports/experiment_of_sample.mako',
                                   affy_id=affy_id,
                                   experiments=experiments,
                                   message=message)

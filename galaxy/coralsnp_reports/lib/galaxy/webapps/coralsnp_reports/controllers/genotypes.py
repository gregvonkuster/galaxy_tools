import calendar
import logging
from datetime import (
    date,
    datetime,
    timedelta
)

import sqlalchemy as sa
from markupsafe import escape

from galaxy.model import corals
from galaxy.webapps.base.controller import (
    BaseUIController,
    web,
)
from galaxy import util
from galaxy.webapps.reports.controllers.jobs import sorter
from galaxy.webapps.coralsnp_reports.controllers.query import ReportQueryBuilder

log = logging.getLogger(__name__)


class Genotypes(BaseUIController, ReportQueryBuilder):

    @web.expose
    def all(self, trans, **kwd):
        message = escape(util.restore_text(kwd.get('message', '')))
        genotypes = []
        for row in trans.sa_session.query(corals.Genotype):
            genotypes.append((row.id, row.coral_mlg_clonal_id, row.coral_mlg_rep_sample_id, row.genetic_coral_species_call))
        return trans.fill_template('/webapps/coralsnp_reports/genotypes.mako', genotypes=genotypes, message=message)

    @web.expose
    def all_by_sample_upload_date(self, trans, **kwd):
        message = escape(util.restore_text(kwd.get('message', '')))
        num_genotypes = trans.sa_session.query(corals.Genotype).count()
        return trans.fill_template('/webapps/coralsnp_reports/genotypes_by_sample_upload_date.mako', num_genotypes=num_genotypes, message=message)

    @web.expose
    def of_sample(self, trans, **kwd):
        message = escape(util.restore_text(kwd.get('message', '')))
        # If specified_date is not received, we'll default to the current month
        affy_id = kwd.get('affy_id')
        genotype_id = kwd.get('genotype_id')
        q = (
            sa.select(
                corals.Genotype.id,
                corals.Genotype.coral_mlg_clonal_id,
                corals.Genotype.coral_mlg_rep_sample_id,
                corals.Genotype.genetic_coral_species_call
            )
            .select_from(corals.Genotype.table)
            .where(corals.Genotype.table.c.id == genotype_id)
            .order_by(corals.Genotype.table.c.id)
        )
        genotypes = []
        for row in trans.sa_session.execute(q):
            genotypes.append((row.coral_mlg_clonal_id, row.coral_mlg_rep_sample_id, row.genetic_coral_species_call))
        return trans.fill_template('/webapps/coralsnp_reports/genotype_of_sample.mako',
                                   affy_id=affy_id,
                                   genotypes=genotypes,
                                   message=message)

    @web.expose
    def per_month(self, trans, **kwd):
        message = escape(util.restore_text(kwd.get('message', '')))
        specs = sorter('date', kwd)
        sort_id = specs.sort_id
        order = specs.order
        arrow = specs.arrow

        q = (
            sa.select(
                self.select_month(corals.Genotype.table.c.create_time).label('date'),
                sa.func.count(corals.Genotype.table.c.id).label('num_genotypes')
            )
            .select_from(corals.Genotype.table)
            .group_by(self.group_by_month(corals.Genotype.table.c.create_time))
            .order_by(specs.exc_order)
        )
        genotypes = []
        for row in trans.sa_session.execute(q):
            genotypes.append((row.date.strftime("%Y-%m"),
                              row.num_genotypes,
                              row.date.strftime("%B"),
                              row.date.strftime("%Y")))
        return trans.fill_template('/webapps/coralsnp_reports/genotypes_per_month.mako',
                                   order=order,
                                   arrow=arrow,
                                   sort_id=sort_id,
                                   genotypes=genotypes,
                                   message=message)

    @web.expose
    def specified_month(self, trans, **kwd):
        message = escape(util.restore_text(kwd.get('message', '')))
        # If specified_date is not received,
        # we'll default to the current month.
        specified_date = kwd.get('specified_date', datetime.utcnow().strftime("%Y-%m-%d"))
        specified_month = specified_date[:7]
        year, month = map(int, specified_month.split("-"))
        start_date = date(year, month, 1)
        end_date = start_date + timedelta(days=calendar.monthrange(year, month)[1])
        month_label = start_date.strftime("%B")
        year_label = start_date.strftime("%Y")
        q = (
            sa.select(
                self.select_day(corals.Genotype.table.c.create_time).label('date'),
                sa.func.count(corals.Genotype.table.c.id).label('num_genotypes')
            )
            .select_from(corals.Genotype.table)
            .where(sa.and_(corals.Genotype.table.c.create_time >= start_date,
                           corals.Genotype.table.c.create_time < end_date))
            .group_by(self.group_by_day(corals.Genotype.table.c.create_time))
            .order_by(sa.desc('date'))
        )
        genotypes = []
        for row in trans.sa_session.execute(q):
            genotypes.append((row.date.strftime("%Y-%m-%d"),
                              row.date.strftime("%d"),
                              row.num_genotypes,
                              row.date.strftime("%A")))
        return trans.fill_template('/webapps/coralsnp_reports/genotypes_specified_month.mako',
                                   month_label=month_label,
                                   year_label=year_label,
                                   month=month,
                                   genotypes=genotypes,
                                   message=message)

    @web.expose
    def specified_date(self, trans, **kwd):
        message = escape(util.restore_text(kwd.get('message', '')))
        # If specified_date is not received, we'll default to the current month
        specified_date = kwd.get('specified_date', datetime.utcnow().strftime("%Y-%m-%d"))
        year, month, day = map(int, specified_date.split("-"))
        start_date = date(year, month, day)
        end_date = start_date + timedelta(days=1)
        day_of_month = start_date.strftime("%d")
        day_label = start_date.strftime("%A")
        month_label = start_date.strftime("%B")
        year_label = start_date.strftime("%Y")
        q = (
            sa.select(
                self.select_day(corals.Genotype.table.c.create_time).label('date'),
                corals.Genotype.id,
                corals.Genotype.coral_mlg_clonal_id,
                corals.Genotype.coral_mlg_rep_sample_id,
                corals.Genotype.genetic_coral_species_call
            )
            .where(sa.and_(corals.Genotype.table.c.create_time >= start_date,
                           corals.Genotype.table.c.create_time < end_date))
            .select_from(corals.Genotype.table)
            .order_by(corals.Genotype.table.c.id)
        )
        genotypes = []
        for row in trans.sa_session.execute(q):
            genotypes.append((row.id, row.coral_mlg_clonal_id, row.coral_mlg_rep_sample_id, row.genetic_coral_species_call))
        return trans.fill_template('/webapps/coralsnp_reports/genotypes_specified_date.mako',
                                   specified_date=start_date,
                                   day_label=day_label,
                                   month_label=month_label,
                                   year_label=year_label,
                                   day_of_month=day_of_month,
                                   genotypes=genotypes,
                                   message=message)

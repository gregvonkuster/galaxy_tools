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

    @web.expose
    def per_month(self, trans, **kwd):
        message = escape(util.restore_text(kwd.get('message', '')))
        specs = sorter('date', kwd)
        sort_id = specs.sort_id
        order = specs.order
        arrow = specs.arrow
        _order = specs.exc_order

        q = sa.select((self.select_month(galaxy.model.corals.Genotype.table.c.create_time).label('date'),
                       sa.func.count(galaxy.model.corals.Genotype.table.c.id).label('num_genotypes')),
                      from_obj=[galaxy.model.corals.Genotype.table],
                      group_by=self.group_by_month(galaxy.model.corals.Genotype.table.c.create_time),
                      order_by=[_order])
        genotypes = []
        for row in q.execute():
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
        # If specified_date is not received, we'll default to the current month
        specified_date = kwd.get('specified_date', datetime.utcnow().strftime("%Y-%m-%d"))
        specified_month = specified_date[:7]
        year, month = map(int, specified_month.split("-"))
        start_date = date(year, month, 1)
        end_date = start_date + timedelta(days=calendar.monthrange(year, month)[1])
        month_label = start_date.strftime("%B")
        year_label = start_date.strftime("%Y")
        q = sa.select((self.select_day(galaxy.model.corals.Genotype.table.c.create_time).label('date'),
                       sa.func.count(galaxy.model.corals.Genotype.table.c.id).label('num_genotypes')),
                      whereclause=sa.and_(galaxy.model.corals.Genotype.table.c.create_time >= start_date,
                                          galaxy.model.corals.Genotype.table.c.create_time < end_date),
                      from_obj=[galaxy.model.corals.Genotype.table],
                      group_by=self.group_by_day(galaxy.model.corals.Genotype.table.c.create_time),
                      order_by=[sa.desc('date')])
        genotypes = []
        for row in q.execute():
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
        q = sa.select((self.select_day(galaxy.model.corals.Genotype.table.c.create_time).label('date'),
                       galaxy.model.corals.Genotype.table.c.id,
                       galaxy.model.corals.Genotype.table.c.coral_mlg_clonal_id,
                       galaxy.model.corals.Genotype.table.c.coral_mlg_rep_sample_id,
                       galaxy.model.corals.Genotype.table.c.genetic_coral_species_call,
                       galaxy.model.corals.Genotype.table.c.bcoral_genet_id),
                      whereclause=sa.and_(galaxy.model.corals.Genotype.table.c.create_time >= start_date,
                                          galaxy.model.corals.Genotype.table.c.create_time < end_date),
                      from_obj=[galaxy.model.corals.Genotype.table],
                      order_by=[galaxy.model.corals.Genotype.table.c.id])
        genotypes = []
        for row in q.execute():
            cols_tup = (row.id, row.coral_mlg_clonal_id, row.coral_mlg_rep_sample_id,
                        row.genetic_coral_species_call, row.bcoral_genet_id)
            genotypes.append((cols_tup))
        return trans.fill_template('/webapps/coralsnp_reports/genotypes_specified_date.mako',
                                   specified_date=start_date,
                                   day_label=day_label,
                                   month_label=month_label,
                                   year_label=year_label,
                                   day_of_month=day_of_month,
                                   genotypes=genotypes,
                                   message=message)

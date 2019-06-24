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


class Samples(BaseUIController, ReportQueryBuilder):

    @web.expose
    def all(self, trans, **kwd):
        message = escape(util.restore_text(kwd.get('message', '')))
        q = sa.select((galaxy.model.corals.Sample.table.c.id,
                       galaxy.model.corals.Sample.table.c.affy_id,
                       galaxy.model.corals.Sample.table.c.sample_id,
                       galaxy.model.corals.Sample.table.c.genotype_id,
                       galaxy.model.corals.Sample.table.c.phenotype_id,
                       galaxy.model.corals.Sample.table.c.experiment_id,
                       galaxy.model.corals.Sample.table.c.colony_id,
                       galaxy.model.corals.Sample.table.c.colony_location,
                       galaxy.model.corals.Sample.table.c.taxonomy_id,
                       galaxy.model.corals.Sample.table.c.collector_id,
                       galaxy.model.corals.Sample.table.c.field_call,
                       galaxy.model.corals.Sample.table.c.collection_date,
                       galaxy.model.corals.Sample.table.c.user_specimen_id,
                       galaxy.model.corals.Sample.table.c.registry_id,
                       galaxy.model.corals.Sample.table.c.depth,
                       galaxy.model.corals.Sample.table.c.dna_extraction_method,
                       galaxy.model.corals.Sample.table.c.dna_concentration,
                       galaxy.model.corals.Sample.table.c.public_after_date,
                       galaxy.model.corals.Sample.table.c.percent_missing_data_coral,
                       galaxy.model.corals.Sample.table.c.percent_reference_coral,
                       galaxy.model.corals.Sample.table.c.percent_alternative_coral,
                       galaxy.model.corals.Sample.table.c.percent_heterozygous_coral),
                      from_obj=[galaxy.model.corals.Sample.table],
                      order_by=[galaxy.model.corals.Sample.table.c.id])
        samples = []
        for row in q.execute():
            try:
                collection_date = row.collection_date.strftime("%Y-%m-%d")
            except Exception:
                collection_date = row.collection_date
            try:
                public_after_date = row.public_after_date.strftime("%Y-%m-%d")
            except Exception:
                public_after_date = row.public_after_date
            cols_tup = (row.affy_id, row.sample_id, row.field_call, row.colony_location,
                        collection_date, row.user_specimen_id, row.registry_id, row.depth,
                        row.dna_extraction_method, row.dna_concentration, public_after_date,
                        row.percent_missing_data_coral, row.percent_reference_coral,
                        row.percent_alternative_coral, row.percent_heterozygous_coral,
                        row.genotype_id, row.phenotype_id, row.collector_id,
                        row.experiment_id, row.colony_id, row.taxonomy_id)
            samples.append(cols_tup)
        return trans.fill_template('/webapps/coralsnp_reports/samples.mako', samples=samples, message=message)

    @web.expose
    def all_by_upload_date(self, trans, **kwd):
        message = escape(util.restore_text(kwd.get('message', '')))
        num_samples = trans.sa_session.query(galaxy.model.corals.Sample).count()
        return trans.fill_template('/webapps/coralsnp_reports/samples_by_upload_date.mako', num_samples=num_samples, message=message)

    @web.expose
    def collected_by(self, trans, **kwd):
        message = escape(util.restore_text(kwd.get('message', '')))
        last_name = kwd.get('last_name')
        first_name = kwd.get('first_name')
        collector_id = kwd.get('collector_id')
        q = sa.select((galaxy.model.corals.Sample.table.c.id,
                       galaxy.model.corals.Sample.table.c.affy_id,
                       galaxy.model.corals.Sample.table.c.sample_id,
                       galaxy.model.corals.Sample.table.c.genotype_id,
                       galaxy.model.corals.Sample.table.c.phenotype_id,
                       galaxy.model.corals.Sample.table.c.experiment_id,
                       galaxy.model.corals.Sample.table.c.colony_id,
                       galaxy.model.corals.Sample.table.c.colony_location,
                       galaxy.model.corals.Sample.table.c.taxonomy_id,
                       galaxy.model.corals.Sample.table.c.field_call,
                       galaxy.model.corals.Sample.table.c.collection_date,
                       galaxy.model.corals.Sample.table.c.user_specimen_id,
                       galaxy.model.corals.Sample.table.c.registry_id,
                       galaxy.model.corals.Sample.table.c.depth,
                       galaxy.model.corals.Sample.table.c.dna_extraction_method,
                       galaxy.model.corals.Sample.table.c.dna_concentration,
                       galaxy.model.corals.Sample.table.c.public_after_date,
                       galaxy.model.corals.Sample.table.c.percent_missing_data_coral,
                       galaxy.model.corals.Sample.table.c.percent_reference_coral,
                       galaxy.model.corals.Sample.table.c.percent_alternative_coral,
                       galaxy.model.corals.Sample.table.c.percent_heterozygous_coral),
                      whereclause=galaxy.model.corals.Sample.table.c.collector_id == collector_id,
                      from_obj=[galaxy.model.corals.Sample.table],
                      order_by=[galaxy.model.corals.Sample.table.c.id])
        samples = []
        for row in q.execute():
            try:
                collection_date = row.collection_date.strftime("%Y-%m-%d")
            except Exception:
                collection_date = row.collection_date
            try:
                public_after_date = row.public_after_date.strftime("%Y-%m-%d")
            except Exception:
                public_after_date = row.public_after_date
            cols_tup = (row.affy_id, row.sample_id, row.field_call, row.colony_location,
                        collection_date, row.user_specimen_id, row.registry_id,
                        row.depth, row.dna_extraction_method, row.dna_concentration,
                        public_after_date, row.percent_missing_data_coral,
                        row.percent_reference_coral, row.percent_alternative_coral,
                        row.percent_heterozygous_coral, row.genotype_id, row.phenotype_id,
                        row.experiment_id, row.colony_id, row.taxonomy_id)
            samples.append(cols_tup)
        return trans.fill_template('/webapps/coralsnp_reports/samples_collected_by.mako',
                                   last_name=last_name,
                                   first_name=first_name,
                                   samples=samples,
                                   message=message)

    @web.expose
    def of_colony(self, trans, **kwd):
        message = escape(util.restore_text(kwd.get('message', '')))
        colony_id = kwd.get('colony_id')
        latitude = kwd.get('latitude')
        longitude = kwd.get('longitude')
        depth = kwd.get('depth')
        reef_id = kwd.get('reef_id')
        q = sa.select((galaxy.model.corals.Sample.table.c.id,
                       galaxy.model.corals.Sample.table.c.affy_id,
                       galaxy.model.corals.Sample.table.c.sample_id,
                       galaxy.model.corals.Sample.table.c.genotype_id,
                       galaxy.model.corals.Sample.table.c.phenotype_id,
                       galaxy.model.corals.Sample.table.c.experiment_id,
                       galaxy.model.corals.Sample.table.c.colony_id,
                       galaxy.model.corals.Sample.table.c.colony_location,
                       galaxy.model.corals.Sample.table.c.taxonomy_id,
                       galaxy.model.corals.Sample.table.c.collector_id,
                       galaxy.model.corals.Sample.table.c.collection_date,
                       galaxy.model.corals.Sample.table.c.user_specimen_id,
                       galaxy.model.corals.Sample.table.c.registry_id,
                       galaxy.model.corals.Sample.table.c.depth,
                       galaxy.model.corals.Sample.table.c.dna_extraction_method,
                       galaxy.model.corals.Sample.table.c.dna_concentration,
                       galaxy.model.corals.Sample.table.c.public_after_date,
                       galaxy.model.corals.Sample.table.c.percent_missing_data_coral,
                       galaxy.model.corals.Sample.table.c.percent_reference_coral,
                       galaxy.model.corals.Sample.table.c.percent_alternative_coral,
                       galaxy.model.corals.Sample.table.c.percent_heterozygous_coral,
                       galaxy.model.corals.Sample.table.c.field_call),
                      whereclause=galaxy.model.corals.Sample.table.c.colony_id == colony_id,
                      from_obj=[galaxy.model.corals.Sample.table],
                      order_by=[galaxy.model.corals.Sample.table.c.id])
        samples = []
        for row in q.execute():
            try:
                collection_date = row.collection_date.strftime("%Y-%m-%d")
            except Exception:
                collection_date = row.collection_date
            try:
                public_after_date = row.public_after_date.strftime("%Y-%m-%d")
            except Exception:
                public_after_date = row.public_after_date
            cols_tup = (row.affy_id, row.sample_id, row.field_call, row.colony_location,
                        collection_date, row.user_specimen_id, row.registry_id, row.depth,
                        row.dna_extraction_method, row.dna_concentration, public_after_date,
                        row.percent_missing_data_coral, row.percent_reference_coral,
                        row.percent_alternative_coral, row.percent_heterozygous_coral,
                        row.genotype_id, row.phenotype_id, row.collector_id, row.experiment_id,
                        row.taxonomy_id)
            samples.append(cols_tup)
        return trans.fill_template('/webapps/coralsnp_reports/samples_of_colony.mako',
                                   latitude=latitude,
                                   longitude=longitude,
                                   depth=depth,
                                   reef_id=reef_id,
                                   samples=samples,
                                   message=message)

    @web.expose
    def of_experiment(self, trans, **kwd):
        message = escape(util.restore_text(kwd.get('message', '')))
        experiment_id = kwd.get('experiment_id')
        seq_facility = kwd.get('seq_facility')
        array_version = kwd.get('array_version')
        result_folder_name = kwd.get('result_folder_name')
        plate_barcode = kwd.get('plate_barcode')
        q = sa.select((galaxy.model.corals.Sample.table.c.id,
                       galaxy.model.corals.Sample.table.c.affy_id,
                       galaxy.model.corals.Sample.table.c.sample_id,
                       galaxy.model.corals.Sample.table.c.genotype_id,
                       galaxy.model.corals.Sample.table.c.phenotype_id,
                       galaxy.model.corals.Sample.table.c.experiment_id,
                       galaxy.model.corals.Sample.table.c.colony_id,
                       galaxy.model.corals.Sample.table.c.colony_location,
                       galaxy.model.corals.Sample.table.c.taxonomy_id,
                       galaxy.model.corals.Sample.table.c.collector_id,
                       galaxy.model.corals.Sample.table.c.collection_date,
                       galaxy.model.corals.Sample.table.c.user_specimen_id,
                       galaxy.model.corals.Sample.table.c.registry_id,
                       galaxy.model.corals.Sample.table.c.depth,
                       galaxy.model.corals.Sample.table.c.dna_extraction_method,
                       galaxy.model.corals.Sample.table.c.dna_concentration,
                       galaxy.model.corals.Sample.table.c.public_after_date,
                       galaxy.model.corals.Sample.table.c.percent_missing_data_coral,
                       galaxy.model.corals.Sample.table.c.percent_reference_coral,
                       galaxy.model.corals.Sample.table.c.percent_alternative_coral,
                       galaxy.model.corals.Sample.table.c.percent_heterozygous_coral,
                       galaxy.model.corals.Sample.table.c.field_call),
                      whereclause=galaxy.model.corals.Sample.table.c.experiment_id == experiment_id,
                      from_obj=[galaxy.model.corals.Sample.table],
                      order_by=[galaxy.model.corals.Sample.table.c.id])
        samples = []
        for row in q.execute():
            try:
                collection_date = row.collection_date.strftime("%Y-%m-%d")
            except Exception:
                collection_date = row.collection_date
            try:
                public_after_date = row.public_after_date.strftime("%Y-%m-%d")
            except Exception:
                public_after_date = row.public_after_date
            cols_tup = (row.affy_id, row.sample_id, row.field_call, row.colony_location,
                        collection_date, row.user_specimen_id, row.registry_id, row.depth,
                        row.dna_extraction_method, row.dna_concentration, public_after_date,
                        row.percent_missing_data_coral, row.percent_reference_coral,
                        row.percent_alternative_coral, row.percent_heterozygous_coral,
                        row.genotype_id, row.phenotype_id, row.collector_id, row.colony_id,
                        row.taxonomy_id)
            samples.append(cols_tup)
        return trans.fill_template('/webapps/coralsnp_reports/samples_of_experiment.mako',
                                   seq_facility=seq_facility,
                                   array_version=array_version,
                                   result_folder_name=result_folder_name,
                                   plate_barcode=plate_barcode,
                                   samples=samples,
                                   message=message)

    @web.expose
    def of_reef(self, trans, **kwd):
        message = escape(util.restore_text(kwd.get('message', '')))
        reef_id = kwd.get('reef_id')
        name = kwd.get('name')
        region = kwd.get('region')
        latitude = kwd.get('latitude')
        longitude = kwd.get('longitude')
        geographic_origin = kwd.get('geographic_origin')
        q = sa.select((galaxy.model.corals.Sample.table.c.id,
                       galaxy.model.corals.Sample.table.c.affy_id,
                       galaxy.model.corals.Sample.table.c.sample_id,
                       galaxy.model.corals.Sample.table.c.genotype_id,
                       galaxy.model.corals.Sample.table.c.phenotype_id,
                       galaxy.model.corals.Sample.table.c.experiment_id,
                       galaxy.model.corals.Sample.table.c.colony_id,
                       galaxy.model.corals.Sample.table.c.colony_location,
                       galaxy.model.corals.Sample.table.c.taxonomy_id,
                       galaxy.model.corals.Sample.table.c.collector_id,
                       galaxy.model.corals.Sample.table.c.collection_date,
                       galaxy.model.corals.Sample.table.c.user_specimen_id,
                       galaxy.model.corals.Sample.table.c.registry_id,
                       galaxy.model.corals.Sample.table.c.depth,
                       galaxy.model.corals.Sample.table.c.dna_extraction_method,
                       galaxy.model.corals.Sample.table.c.dna_concentration,
                       galaxy.model.corals.Sample.table.c.public_after_date,
                       galaxy.model.corals.Sample.table.c.percent_missing_data_coral,
                       galaxy.model.corals.Sample.table.c.percent_reference_coral,
                       galaxy.model.corals.Sample.table.c.percent_alternative_coral,
                       galaxy.model.corals.Sample.table.c.percent_heterozygous_coral,
                       galaxy.model.corals.Sample.table.c.field_call,
                       galaxy.model.corals.Colony.table.c.reef_id,
                       galaxy.model.corals.Reef.table.c.id),
                      from_obj=[galaxy.model.corals.Sample.table,
                                galaxy.model.corals.Colony.table,
                                galaxy.model.corals.Reef.table],
                      whereclause=sa.and_(galaxy.model.corals.Sample.table.c.colony_id == galaxy.model.corals.Colony.table.c.id,
                                          galaxy.model.corals.Colony.table.c.reef_id == galaxy.model.corals.Reef.table.c.id,
                                          galaxy.model.corals.Reef.table.c.id == reef_id), 
                      order_by=[galaxy.model.corals.Sample.table.c.id])
        samples = []
        for row in q.execute():
            try:
                collection_date = row.collection_date.strftime("%Y-%m-%d")
            except Exception:
                collection_date = row.collection_date
            try:
                public_after_date = row.public_after_date.strftime("%Y-%m-%d")
            except Exception:
                public_after_date = row.public_after_date
            cols_tup = (row.affy_id, row.sample_id, row.field_call, row.colony_location,
                        collection_date, row.user_specimen_id, row.registry_id, row.depth,
                        row.dna_extraction_method, row.dna_concentration, public_after_date,
                        row.percent_missing_data_coral, row.percent_reference_coral,
                        row.percent_alternative_coral, row.percent_heterozygous_coral,
                        row.genotype_id, row.phenotype_id, row.collector_id, row.experiment_id,
                        row.colony_id, row.taxonomy_id)
            samples.append(cols_tup)
        return trans.fill_template('/webapps/coralsnp_reports/samples_of_reef.mako',
                                   name=name,
                                   region=region,
                                   latitude=latitude,
                                   longitude=longitude,
                                   geographic_origin=geographic_origin,
                                   samples=samples,
                                   message=message)

    @web.expose
    def of_taxonomy(self, trans, **kwd):
        message = escape(util.restore_text(kwd.get('message', '')))
        taxonomy_id = kwd.get('taxonomy_id')
        species_name = kwd.get('species_name')
        genus_name = kwd.get('genus_name')
        q = sa.select((galaxy.model.corals.Sample.table.c.id,
                       galaxy.model.corals.Sample.table.c.affy_id,
                       galaxy.model.corals.Sample.table.c.sample_id,
                       galaxy.model.corals.Sample.table.c.genotype_id,
                       galaxy.model.corals.Sample.table.c.phenotype_id,
                       galaxy.model.corals.Sample.table.c.experiment_id,
                       galaxy.model.corals.Sample.table.c.colony_id,
                       galaxy.model.corals.Sample.table.c.colony_location,
                       galaxy.model.corals.Sample.table.c.taxonomy_id,
                       galaxy.model.corals.Sample.table.c.collector_id,
                       galaxy.model.corals.Sample.table.c.collection_date,
                       galaxy.model.corals.Sample.table.c.user_specimen_id,
                       galaxy.model.corals.Sample.table.c.registry_id,
                       galaxy.model.corals.Sample.table.c.depth,
                       galaxy.model.corals.Sample.table.c.dna_extraction_method,
                       galaxy.model.corals.Sample.table.c.dna_concentration,
                       galaxy.model.corals.Sample.table.c.public_after_date,
                       galaxy.model.corals.Sample.table.c.percent_missing_data_coral,
                       galaxy.model.corals.Sample.table.c.percent_reference_coral,
                       galaxy.model.corals.Sample.table.c.percent_alternative_coral,
                       galaxy.model.corals.Sample.table.c.percent_heterozygous_coral,
                       galaxy.model.corals.Sample.table.c.field_call),
                      whereclause=galaxy.model.corals.Sample.table.c.taxonomy_id == taxonomy_id,
                      from_obj=[galaxy.model.corals.Sample.table],
                      order_by=[galaxy.model.corals.Sample.table.c.id])
        samples = []
        for row in q.execute():
            try:
                collection_date = row.collection_date.strftime("%Y-%m-%d")
            except Exception:
                collection_date = row.collection_date
            try:
                public_after_date = row.public_after_date.strftime("%Y-%m-%d")
            except Exception:
                public_after_date = row.public_after_date
            cols_tup = (row.affy_id, row.sample_id, row.field_call, row.colony_location,
                        collection_date, row.user_specimen_id, row.registry_id, row.depth,
                        row.dna_extraction_method, row.dna_concentration, public_after_date,
                        row.percent_missing_data_coral, row.percent_reference_coral,
                        row.percent_alternative_coral, row.percent_heterozygous_coral,
                        row.genotype_id, row.phenotype_id, row.collector_id, row.experiment_id,
                        row.colony_id)
            samples.append(cols_tup)
        return trans.fill_template('/webapps/coralsnp_reports/samples_of_taxonomy.mako',
                                   species_name=species_name,
                                   genus_name=genus_name,
                                   samples=samples,
                                   message=message)

    @web.expose
    def per_month(self, trans, **kwd):
        message = escape(util.restore_text(kwd.get('message', '')))
        specs = sorter('date', kwd)
        sort_id = specs.sort_id
        order = specs.order
        arrow = specs.arrow
        _order = specs.exc_order

        q = sa.select((self.select_month(galaxy.model.corals.Sample.table.c.create_time).label('date'),
                       sa.func.count(galaxy.model.corals.Sample.table.c.id).label('num_samples')),
                      from_obj=[galaxy.model.corals.Sample.table],
                      group_by=self.group_by_month(galaxy.model.corals.Sample.table.c.create_time),
                      order_by=[_order])
        samples = []
        for row in q.execute():
            samples.append((row.date.strftime("%Y-%m"),
                            row.num_samples,
                            row.date.strftime("%B"),
                            row.date.strftime("%Y")))
        return trans.fill_template('/webapps/coralsnp_reports/samples_per_month.mako',
                                   order=order,
                                   arrow=arrow,
                                   sort_id=sort_id,
                                   samples=samples,
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
        q = sa.select((self.select_day(galaxy.model.corals.Sample.table.c.create_time).label('date'),
                       galaxy.model.corals.Sample.table.c.id,
                       galaxy.model.corals.Sample.table.c.affy_id,
                       galaxy.model.corals.Sample.table.c.sample_id,
                       galaxy.model.corals.Sample.table.c.genotype_id,
                       galaxy.model.corals.Sample.table.c.phenotype_id,
                       galaxy.model.corals.Sample.table.c.experiment_id,
                       galaxy.model.corals.Sample.table.c.colony_id,
                       galaxy.model.corals.Sample.table.c.colony_location,
                       galaxy.model.corals.Sample.table.c.taxonomy_id,
                       galaxy.model.corals.Sample.table.c.collector_id,
                       galaxy.model.corals.Sample.table.c.field_call,
                       galaxy.model.corals.Sample.table.c.collection_date,
                       galaxy.model.corals.Sample.table.c.user_specimen_id,
                       galaxy.model.corals.Sample.table.c.registry_id,
                       galaxy.model.corals.Sample.table.c.depth,
                       galaxy.model.corals.Sample.table.c.dna_extraction_method,
                       galaxy.model.corals.Sample.table.c.dna_concentration,
                       galaxy.model.corals.Sample.table.c.public_after_date,
                       galaxy.model.corals.Sample.table.c.percent_missing_data_coral,
                       galaxy.model.corals.Sample.table.c.percent_reference_coral,
                       galaxy.model.corals.Sample.table.c.percent_alternative_coral,
                       galaxy.model.corals.Sample.table.c.percent_heterozygous_coral),
                      whereclause=sa.and_(galaxy.model.corals.Sample.table.c.create_time >= start_date,
                                          galaxy.model.corals.Sample.table.c.create_time < end_date),
                      from_obj=[galaxy.model.corals.Sample.table],
                      order_by=[galaxy.model.corals.Sample.table.c.id])
        samples = []
        for row in q.execute():
            try:
                collection_date = row.collection_date.strftime("%Y-%m-%d")
            except Exception:
                collection_date = row.collection_date
            try:
                public_after_date = row.public_after_date.strftime("%Y-%m-%d")
            except Exception:
                public_after_date = row.public_after_date
            cols_tup = (row.affy_id, row.sample_id, row.field_call, row.colony_location,
                        collection_date, row.user_specimen_id,
                        row.registry_id, row.depth, row.dna_extraction_method,
                        row.dna_concentration, public_after_date, row.percent_missing_data_coral,
                        row.percent_reference_coral, row.percent_alternative_coral,
                        row.percent_heterozygous_coral, row.genotype_id, row.phenotype_id,
                        row.collector_id, row.experiment_id, row.colony_id, row.taxonomy_id)
            samples.append(cols_tup)
        return trans.fill_template('/webapps/coralsnp_reports/samples_specified_date.mako',
                                   specified_date=start_date,
                                   day_label=day_label,
                                   month_label=month_label,
                                   year_label=year_label,
                                   day_of_month=day_of_month,
                                   samples=samples,
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
        q = sa.select((self.select_day(galaxy.model.corals.Sample.table.c.create_time).label('date'),
                       sa.func.count(galaxy.model.corals.Sample.table.c.id).label('num_samples')),
                      whereclause=sa.and_(galaxy.model.corals.Sample.table.c.create_time >= start_date,
                                          galaxy.model.corals.Sample.table.c.create_time < end_date),
                      from_obj=[galaxy.model.corals.Sample.table],
                      group_by=self.group_by_day(galaxy.model.corals.Sample.table.c.create_time),
                      order_by=[sa.desc('date')])
        samples = []
        for row in q.execute():
            samples.append((row.date.strftime("%Y-%m-%d"),
                            row.date.strftime("%d"),
                            row.num_samples,
                            row.date.strftime("%A")))
        return trans.fill_template('/webapps/coralsnp_reports/samples_specified_month.mako',
                                   month_label=month_label,
                                   year_label=year_label,
                                   month=month,
                                   samples=samples,
                                   message=message)

    @web.expose
    def with_genotype(self, trans, **kwd):
        message = escape(util.restore_text(kwd.get('message', '')))
        genotype_id = kwd.get('genotype_id')
        coral_mlg_clonal_id = kwd.get('coral_mlg_clonal_id')
        coral_mlg_rep_sample_id = kwd.get('coral_mlg_rep_sample_id')
        genetic_coral_species_call = kwd.get('genetic_coral_species_call')
        bcoral_genet_id = kwd.get('bcoral_genet_id')
        q = sa.select((galaxy.model.corals.Sample.table.c.id,
                       galaxy.model.corals.Sample.table.c.affy_id,
                       galaxy.model.corals.Sample.table.c.sample_id,
                       galaxy.model.corals.Sample.table.c.phenotype_id,
                       galaxy.model.corals.Sample.table.c.experiment_id,
                       galaxy.model.corals.Sample.table.c.colony_id,
                       galaxy.model.corals.Sample.table.c.colony_location,
                       galaxy.model.corals.Sample.table.c.taxonomy_id,
                       galaxy.model.corals.Sample.table.c.collector_id,
                       galaxy.model.corals.Sample.table.c.collection_date,
                       galaxy.model.corals.Sample.table.c.user_specimen_id,
                       galaxy.model.corals.Sample.table.c.registry_id,
                       galaxy.model.corals.Sample.table.c.depth,
                       galaxy.model.corals.Sample.table.c.dna_extraction_method,
                       galaxy.model.corals.Sample.table.c.dna_concentration,
                       galaxy.model.corals.Sample.table.c.public_after_date,
                       galaxy.model.corals.Sample.table.c.percent_missing_data_coral,
                       galaxy.model.corals.Sample.table.c.percent_reference_coral,
                       galaxy.model.corals.Sample.table.c.percent_alternative_coral,
                       galaxy.model.corals.Sample.table.c.percent_heterozygous_coral,
                       galaxy.model.corals.Sample.table.c.field_call),
                      whereclause=galaxy.model.corals.Sample.table.c.genotype_id == genotype_id,
                      from_obj=[galaxy.model.corals.Sample.table],
                      order_by=[galaxy.model.corals.Sample.table.c.id])
        samples = []
        for row in q.execute():
            try:
                collection_date = row.collection_date.strftime("%Y-%m-%d")
            except Exception:
                collection_date = row.collection_date
            try:
                public_after_date = row.public_after_date.strftime("%Y-%m-%d")
            except Exception:
                public_after_date = row.public_after_date
            cols_tup = (row.affy_id, row.sample_id, row.field_call, row.colony_location,
                        collection_date, row.user_specimen_id, row.registry_id, row.depth,
                        row.dna_extraction_method, row.dna_concentration, public_after_date,
                        row.percent_missing_data_coral, row.percent_reference_coral,
                        row.percent_alternative_coral, row.percent_heterozygous_coral,
                        row.phenotype_id, row.collector_id, row.experiment_id, row.colony_id,
                        row.taxonomy_id)
            samples.append(cols_tup)
        return trans.fill_template('/webapps/coralsnp_reports/samples_with_genotype.mako',
                                   coral_mlg_clonal_id=coral_mlg_clonal_id,
                                   coral_mlg_rep_sample_id=coral_mlg_rep_sample_id,
                                   genetic_coral_species_call=genetic_coral_species_call,
                                   bcoral_genet_id=bcoral_genet_id,
                                   samples=samples,
                                   message=message)

    @web.expose
    def with_phenotype(self, trans, **kwd):
        message = escape(util.restore_text(kwd.get('message', '')))
        phenotype_id = kwd.get('phenotype_id')
        disease_resist = kwd.get('disease_resist')
        bleach_resist = kwd.get('bleach_resist')
        mortality = kwd.get('mortality')
        tle = kwd.get('tle')
        spawning = kwd.get('spawning')
        sperm_motility = kwd.get('sperm_motility')
        healing_time = kwd.get('healing_time')
        q = sa.select((galaxy.model.corals.Sample.table.c.id,
                       galaxy.model.corals.Sample.table.c.affy_id,
                       galaxy.model.corals.Sample.table.c.sample_id,
                       galaxy.model.corals.Sample.table.c.genotype_id,
                       galaxy.model.corals.Sample.table.c.experiment_id,
                       galaxy.model.corals.Sample.table.c.colony_id,
                       galaxy.model.corals.Sample.table.c.colony_location,
                       galaxy.model.corals.Sample.table.c.taxonomy_id,
                       galaxy.model.corals.Sample.table.c.collector_id,
                       galaxy.model.corals.Sample.table.c.collection_date,
                       galaxy.model.corals.Sample.table.c.user_specimen_id,
                       galaxy.model.corals.Sample.table.c.registry_id,
                       galaxy.model.corals.Sample.table.c.depth,
                       galaxy.model.corals.Sample.table.c.dna_extraction_method,
                       galaxy.model.corals.Sample.table.c.dna_concentration,
                       galaxy.model.corals.Sample.table.c.public_after_date,
                       galaxy.model.corals.Sample.table.c.percent_missing_data_coral,
                       galaxy.model.corals.Sample.table.c.percent_reference_coral,
                       galaxy.model.corals.Sample.table.c.percent_alternative_coral,
                       galaxy.model.corals.Sample.table.c.percent_heterozygous_coral,
                       galaxy.model.corals.Sample.table.c.field_call),
                      whereclause=galaxy.model.corals.Sample.table.c.phenotype_id == phenotype_id,
                      from_obj=[galaxy.model.corals.Sample.table],
                      order_by=[galaxy.model.corals.Sample.table.c.id])
        samples = []
        for row in q.execute():
            try:
                collection_date = row.collection_date.strftime("%Y-%m-%d")
            except Exception:
                collection_date = row.collection_date
            try:
                public_after_date = row.public_after_date.strftime("%Y-%m-%d")
            except Exception:
                public_after_date = row.public_after_date
            cols_tup = (row.affy_id, row.sample_id, row.field_call, row.colony_location,
                        collection_date, row.user_specimen_id, row.registry_id, row.depth,
                        row.dna_extraction_method, row.dna_concentration, public_after_date,
                        row.percent_missing_data_coral, row.percent_reference_coral,
                        row.percent_alternative_coral, row.percent_heterozygous_coral,
                        row.genotype_id, row.collector_id, row.experiment_id, row.colony_id,
                        row.taxonomy_id)
            samples.append(cols_tup)
        return trans.fill_template('/webapps/coralsnp_reports/samples_with_phenotype.mako',
                                   disease_resist=disease_resist,
                                   bleach_resist=bleach_resist,
                                   mortality=mortality,
                                   tle=tle,
                                   spawning=spawning,
                                   sperm_motility=sperm_motility,
                                   healing_time=healing_time,
                                   samples=samples,
                                   message=message)

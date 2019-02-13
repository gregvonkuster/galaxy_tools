"""
Details of how the corals data model objects are mapped onto the relational database
are encapsulated here.
"""
import datetime
import logging

from sqlalchemy import (
    Boolean,
    Column,
    DateTime,
    ForeignKey,
    Integer,
    MetaData,
    Numeric,
    String,
    Table
)
from sqlalchemy.orm import mapper, relation

import galaxy.model.corals as corals_model
from galaxy.model.base import ModelMapping
from galaxy.model.custom_types import TrimmedString
from galaxy.model.orm.engine_factory import build_engine
from galaxy.model.orm.now import now

log = logging.getLogger(__name__)

metadata = MetaData()

# Get current date plus one year for insertion into
# the public_after_date column of the sample table.
today = datetime.date.today()
try:
    # Return the same day of the year.
    year_from_now = today.replace(year=today.year+1)
except Exception:
    # Handle leap years.
    year_from_now = today + (datetime.date(today.year + 1, 1, 1) - datetime.date(today.year, 1, 1))

corals_model.Collector.table = Table("collector", metadata,
    Column("id", Integer, primary_key=True),
    Column("create_time", DateTime, default=now),
    Column("update_time", DateTime, default=now, onupdate=now),
    Column("person_id", Integer, ForeignKey("person.id"), index=True),
    Column("contact_id", Integer, ForeignKey("person.id"), index=True))

corals_model.Colony.table = Table("colony", metadata,
    Column("id", Integer, primary_key=True),
    Column("create_time", DateTime, default=now),
    Column("update_time", DateTime, default=now, onupdate=now),
    Column("latitude", Numeric(15, 6)),
    Column("longitude", Numeric(15, 6)),
    Column("depth", Integer),
    Column("reef_id", Integer, ForeignKey("reef.id"), index=True))

corals_model.Experiment.table = Table("experiment", metadata,
    Column("id", Integer, primary_key=True),
    Column("create_time", DateTime, default=now),
    Column("update_time", DateTime, default=now, onupdate=now),
    Column("seq_facility", String),
    Column("array_version", TrimmedString(255)))

corals_model.Fragment.table = Table("fragment", metadata,
    Column("id", Integer, primary_key=True),
    Column("create_time", DateTime, default=now),
    Column("update_time", DateTime, default=now, onupdate=now),
    Column("colony_id", Integer, ForeignKey("colony.id"), index=True))

corals_model.Genotype.table = Table("genotype", metadata,
    Column("id", Integer, primary_key=True),
    Column("create_time", DateTime, default=now),
    Column("update_time", DateTime, default=now, onupdate=now),
    Column("coral_mlg_clonal_id", TrimmedString(255)),
    Column("symbio_mlg_clonal_id", TrimmedString(255)),
    Column("genetic_coral_species_call", TrimmedString(255)))

corals_model.Person.table = Table("person", metadata,
    Column("id", Integer, primary_key=True),
    Column("create_time", DateTime, default=now),
    Column("update_time", DateTime, default=now, onupdate=now),
    Column("last_name", TrimmedString(255)),
    Column("first_name", TrimmedString(255)),
    Column("organization", TrimmedString(255)),
    Column("email", TrimmedString(255)))

corals_model.Probe_annotation.table = Table("probe_annotation", metadata,
    Column("id", Integer, primary_key=True),
    Column("create_time", DateTime, default=now),
    Column("update_time", DateTime, default=now, onupdate=now),
    Column("probe_set_id", TrimmedString(255)),
    Column("affy_snp_id", TrimmedString(255)),
    Column("chr_id", Integer),
    Column("start", Integer),
    Column("strand", TrimmedString(255)),
    Column("flank", TrimmedString(255)),
    Column("allele_a", TrimmedString(255)),
    Column("allele_b", TrimmedString(255)),
    Column("allele_frequencies", TrimmedString(255)),
    Column("annotation_notes", TrimmedString(255)),
    Column("allele_count", TrimmedString(255)),
    Column("ordered_alleles", TrimmedString(255)),
    Column("chrtype", TrimmedString(255)),
    Column("custchr", TrimmedString(255)),
    Column("custid", TrimmedString(255)),
    Column("custpos", TrimmedString(255)),
    Column("organism", TrimmedString(255)),
    Column("pconvert", TrimmedString(255)),
    Column("recommendation", TrimmedString(255)),
    Column("refstr", TrimmedString(255)),
    Column("snppriority", TrimmedString(255)))

corals_model.Reef.table = Table("reef", metadata,
    Column("id", Integer, primary_key=True),
    Column("create_time", DateTime, default=now),
    Column("update_time", DateTime, default=now, onupdate=now),
    Column("name", TrimmedString(255)),
    Column("region", TrimmedString(255)),
    Column("latitude", Numeric(15, 6)),
    Column("longitude", Numeric(15, 6)))

corals_model.Phenotype.table = Table("phenotype", metadata,
    Column("id", Integer, primary_key=True),
    Column("create_time", DateTime, default=now),
    Column("update_time", DateTime, default=now, onupdate=now),
    Column("disease_resist", TrimmedString(255)),
    Column("bleach_resist", TrimmedString(255)),
    Column("mortality", TrimmedString(255)),
    Column("tle", TrimmedString(255)),
    Column("spawning", Boolean),
    Column("sperm_motility", Numeric(15, 6)),
    Column("healing_time", Numeric(15, 6)))

corals_model.Sample.table = Table("sample", metadata,
   Column("id", Integer, primary_key=True),
    Column("create_time", DateTime, default=now),
    Column("update_time", DateTime, default=now, onupdate=now),
    Column("affy_id", TrimmedString(255), index=True, nullable=False),
    Column("sample_id", TrimmedString(255), index=True, nullable=False),
    Column("genotype_id", Integer, ForeignKey("genotype.id"), index=True),
    Column("phenotype_id", Integer, ForeignKey("phenotype.id"), index=True),
    Column("experiment_id", Integer, ForeignKey("experiment.id"), index=True),
    Column("colony_id", Integer, ForeignKey("colony.id"), index=True),
    Column("colony_location", TrimmedString(255)),
    Column("fragment_id", Integer, ForeignKey("fragment.id"), index=True),
    Column("taxonomy_id", Integer, ForeignKey("taxonomy.id"), index=True),
    Column("collector_id", Integer, ForeignKey("collector.id"), index=True),
    Column("collection_date", DateTime),
    Column("user_specimen_id", TrimmedString(255)),
    Column("depth", Integer),
    Column("dna_extraction_method", TrimmedString(255)),
    Column("dna_concentration", Numeric(10, 6)),
    Column("public", Boolean),
    Column("public_after_date", DateTime, default=year_from_now),
    Column("percent_missing_data_coral", Numeric(15, 6)),
    Column("percent_missing_data_sym", Numeric(15, 6)),
    Column("percent_reference_coral", Numeric(15, 6)),
    Column("percent_reference_sym", Numeric(15, 6)),
    Column("percent_alternative_coral", Numeric(15, 6)),
    Column("percent_alternative_sym", Numeric(15, 6)),
    Column("percent_hererozygous_coral", Numeric(15, 6)),
    Column("percent_hererozygous_sym", Numeric(15, 6)))

corals_model.Taxonomy.table = Table("taxonomy", metadata,
    Column("id", Integer, primary_key=True),
    Column("create_time", DateTime, default=now),
    Column("update_time", DateTime, default=now, onupdate=now),
    Column("species_name", TrimmedString(255)),
    Column("genus_name", TrimmedString(255)))

# With the tables defined we can define the mappers and setup the
# relationships between the model objects.
mapper(corals_model.Collector, corals_model.Collector.table, properties=dict(
    person=relation(corals_model.Person,
                    primaryjoin=(corals_model.Collector.table.c.person_id == corals_model.Person.table.c.id)),
    contact=relation(corals_model.Person,
                    primaryjoin=(corals_model.Collector.table.c.contact_id == corals_model.Person.table.c.id))))

mapper(corals_model.Colony, corals_model.Colony.table, properties=dict(
    reef=relation(corals_model.Reef,
                  lazy=False,
                  backref="colonies",
                  primaryjoin=(corals_model.Colony.table.c.reef_id == corals_model.Reef.table.c.id))))

mapper(corals_model.Experiment, corals_model.Experiment.table, properties=None)

mapper(corals_model.Fragment, corals_model.Fragment.table, properties=dict(
    colony=relation(corals_model.Colony,
                    lazy=False,
                    backref="fragments",
                    primaryjoin=(corals_model.Fragment.table.c.colony_id == corals_model.Colony.table.c.id))))

mapper(corals_model.Genotype, corals_model.Genotype.table, properties=None)

mapper(corals_model.Person, corals_model.Person.table, properties=None)

mapper(corals_model.Phenotype, corals_model.Phenotype.table, properties=None)

mapper(corals_model.Probe_annotation, corals_model.Probe_annotation.table, properties=None)

mapper(corals_model.Reef, corals_model.Reef.table, properties=None)

mapper(corals_model.Sample, corals_model.Sample.table, properties=dict(
    genotype=relation(corals_model.Genotype,
                      lazy=False,
                      backref="samples",
                      primaryjoin=(corals_model.Sample.table.c.genotype_id == corals_model.Genotype.table.c.id)),
    phenotype=relation(corals_model.Phenotype,
                       lazy=False,
                       backref="samples",
                       primaryjoin=(corals_model.Sample.table.c.phenotype_id == corals_model.Phenotype.table.c.id)),
    experiment=relation(corals_model.Experiment,
                        lazy=False,
                        backref="samples",
                        primaryjoin=(corals_model.Sample.table.c.experiment_id == corals_model.Experiment.table.c.id)),
    fragment=relation(corals_model.Fragment,
                      lazy=False,
                      backref="samples",
                      primaryjoin=(corals_model.Sample.table.c.fragment_id == corals_model.Fragment.table.c.id)),
    colony=relation(corals_model.Colony,
                    lazy=False,
                    backref="matching_samples",
                    primaryjoin=(corals_model.Sample.table.c.colony_id == corals_model.Colony.table.c.id)),
    taxonomy=relation(corals_model.Taxonomy,
                      lazy=False,
                      backref="samples",
                      primaryjoin=(corals_model.Sample.table.c.taxonomy_id == corals_model.Taxonomy.table.c.id)),
    collector=relation(corals_model.Collector,
                       lazy=False,
                       backref="samples",
                       primaryjoin=(corals_model.Sample.table.c.collector_id == corals_model.Collector.table.c.id))))

mapper(corals_model.Taxonomy, corals_model.Taxonomy.table, properties=None)


def init(url, engine_options={}, create_tables=False):
    """Connect mappings to the database"""
    # Load the appropriate db module
    engine = build_engine(url, engine_options)
    # Connect the metadata to the database.
    metadata.bind = engine
    result = ModelMapping([corals_model], engine=engine)
    # Create tables if needed
    if create_tables:
        metadata.create_all()
        # metadata.engine.commit()
    result.create_tables = create_tables
    # load local galaxy security policy
    return result

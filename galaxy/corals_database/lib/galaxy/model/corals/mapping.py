"""
Details of how the corals data model objects are mapped onto the relational database
are encapsulated here.
"""

import logging

from sqlalchemy import (
    and_,
    asc,
    Boolean,
    Column,
    DateTime,
    desc,
    false,
    ForeignKey,
    func,
    Integer,
    MetaData,
    not_,
    Numeric,
    select,
    String, Table,
    TEXT,
    Text,
    true,
    Unicode,
    UniqueConstraint,
    VARCHAR
)
from sqlalchemy.ext.associationproxy import association_proxy
from sqlalchemy.ext.orderinglist import ordering_list
from sqlalchemy.orm import backref, class_mapper, column_property, deferred, mapper, object_session, relation
from sqlalchemy.orm.collections import attribute_mapped_collection
from sqlalchemy.sql import exists
from sqlalchemy.types import BigInteger

import galaxy.model.corals as corals_model
from galaxy.model.base import ModelMapping
from galaxy.model.custom_types import JSONType, MetadataType, TrimmedString, UUIDType
from galaxy.model.orm.engine_factory import build_engine
from galaxy.model.orm.now import now
from galaxy.security import GalaxyRBACAgent

log = logging.getLogger(__name__)

metadata = MetaData()


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
    Column("latitude", Numeric(15, 10)),
    Column("longitude", Numeric(15, 10)),
    Column("depth", Integer),
    Column("reef_id", Integer, ForeignKey("reef.id"), index=True))

corals_model.Coralvcf_allele.table = Table("coralvcf_allele", metadata,
    Column("id", Integer, primary_key=True),
    Column("create_time", DateTime, default=now),
    Column("update_time", DateTime, default=now, onupdate=now),
    Column("chr", TrimmedString(255)),
    Column("pos", Integer),
    Column("sample_id", ForeignKey("sample.id"), index=True),
    Column("ref", TrimmedString(255)),
    Column("alt", TrimmedString(255)),
    Column("qual", Numeric(10, 5)),
    Column("filter", TrimmedString(255)),
    Column("ac", Integer),
    Column("an", Integer),
    Column("bqb", Numeric(10, 5)),
    Column("dp", Integer),
    Column("hob", Numeric(10, 5)),
    Column("icb", Numeric(10, 5)),
    Column("idf", Integer),
    Column("imf", Numeric(10, 5)),
    Column("indel", Boolean),
    Column("mq", Integer),
    Column("mqof", Numeric(10, 5)),
    Column("mqb", Numeric(10, 5)),
    Column("mqsb", Numeric(10, 5)),
    Column("rpb", Numeric(10, 5)),
    Column("sgb", Numeric(10, 5)),
    Column("vdb", Numeric(10, 5)))

corals_model.Experiment.table = Table("experiment", metadata,
    Column("id", Integer, primary_key=True),
    Column("create_time", DateTime, default=now),
    Column("update_time", DateTime, default=now, onupdate=now),
    Column("seq_facility", String),
    Column("array_version", TrimmedString(255)),
    Column("data_sharing", TrimmedString(255)),
    Column("data_hold", TrimmedString(255)))

corals_model.Genotype.table = Table("genotype", metadata,
    Column("id", Integer, primary_key=True),
    Column("create_time", DateTime, default=now),
    Column("update_time", DateTime, default=now, onupdate=now),
    Column("coral_mlg_clonal_id", TrimmedString(255)),
    Column("symbiot_mlg_clonal_id", TrimmedString(255)),
    Column("genetic_coral_species_call", TrimmedString(255)),
    Column("percent_missing_data", Numeric(10, 5)),
    Column("percent_apalm", Numeric(10, 5)),
    Column("percent_acerv", Numeric(10, 5)),
    Column("percent_mixed", Numeric(10, 5)))

corals_model.Idx_annotation.table = Table("idx_annotation", metadata,
    Column("id", Integer, primary_key=True),
    Column("create_time", DateTime, default=now),
    Column("update_time", DateTime, default=now, onupdate=now),
    Column("chip_version", TrimmedString(255)),
    Column("probe_set_id", TrimmedString(255)),
    Column("chr_id", Integer),
    Column("start", Integer),
    Column("stop", Integer),
    Column("strand", TrimmedString(255)),
    Column("dbsnp_rs_id", TrimmedString(255)),
    Column("strand_vs_dbsnp", TrimmedString(255)),
    Column("probe_count", Integer),
    Column("cytoband", TrimmedString(255)),
    Column("chrx_par", Integer),
    Column("allele_a", TrimmedString(255)),
    Column("allele_b", TrimmedString(255)))

corals_model.Person.table = Table("person", metadata,
    Column("id", Integer, primary_key=True),
    Column("create_time", DateTime, default=now),
    Column("update_time", DateTime, default=now, onupdate=now),
    Column("lastname", TrimmedString(255)),
    Column("firstname", TrimmedString(255)),
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
    Column("ref_str", TrimmedString(255)),
    Column("snp_priority", TrimmedString(255)))

corals_model.Reef.table = Table("reef", metadata,
    Column("id", Integer, primary_key=True),
    Column("create_time", DateTime, default=now),
    Column("update_time", DateTime, default=now, onupdate=now),
    Column("name", TrimmedString(255)),
    Column("region", TrimmedString(255)),
    Column("latitude", Numeric(15, 10)),
    Column("longitude", Numeric(15, 10)))

corals_model.Sample.table = Table("sample", metadata,
    Column("id", Integer, primary_key=True),
    Column("create_time", DateTime, default=now),
    Column("update_time", DateTime, default=now, onupdate=now),
    Column("sample_id", TrimmedString(255), index=True, nullable=False),
    Column("genotype_id", Integer, ForeignKey("genotype.id"), index=True),
    Column("experiment_id", Integer, ForeignKey("experiment.id"), index=True),
    Column("colony_id", Integer, ForeignKey("colony.id"), index=True),
    Column("taxonomy_id", Integer, ForeignKey("taxonomy.id"), index=True),
    Column("collector_id", Integer, ForeignKey("collector.id"), index=True),
    Column("collection_date", DateTime),
    Column("user_specimen_id", TrimmedString(255)),
    Column("depth", Integer),
    Column("dna_extraction_method", TrimmedString(255)),
    Column("dna_concentration", Numeric(10, 5)),
    Column("duplicate_sample", Boolean),
    Column("public", Boolean))

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

mapper(corals_model.Coralvcf_allele, corals_model.Coralvcf_allele.table, properties=dict(
    sample=relation(corals_model.Sample,
                  lazy=False,
                  backref="coralvcf_alleles",
                  primaryjoin=(corals_model.Coralvcf_allele.table.c.sample_id == corals_model.Sample.table.c.id))))

mapper(corals_model.Experiment, corals_model.Experiment.table, properties=None)

mapper(corals_model.Genotype, corals_model.Genotype.table, properties=None)

mapper(corals_model.Idx_annotation, corals_model.Idx_annotation.table, properties=None)

mapper(corals_model.Person, corals_model.Person.table, properties=None)

mapper(corals_model.Probe_annotation, corals_model.Probe_annotation.table, properties=None)

mapper(corals_model.Reef, corals_model.Reef.table, properties=None)

mapper(corals_model.Sample, corals_model.Sample.table, properties=dict(
    genotype=relation(corals_model.Genotype,
                      lazy=False,
                      backref="samples",
                      primaryjoin=(corals_model.Sample.table.c.genotype_id == corals_model.Genotype.table.c.id)),
    experiment=relation(corals_model.Experiment,
                        lazy=False,
                        backref="samples",
                        primaryjoin=(corals_model.Sample.table.c.experiment_id == corals_model.Experiment.table.c.id)),
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

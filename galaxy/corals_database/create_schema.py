
from datetime import datetime
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
    String,
    Table,
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
from sqlalchemy.types import BigInteger, TypeDecorator
from sqlalchemy import MetaData
from sqlalchemy import create_engine

now = datetime.utcnow
database_file = "coralSNPs.sqlite"
database_connection = "sqlite:///./%s?isolation_level=IMMEDIATE" % database_file
engine = create_engine(database_connection)
metadata = MetaData(bind=engine)


class TrimmedString(TypeDecorator):
    impl = String

    def process_bind_param(self, value, dialect):
        """Automatically truncate string values"""
        if self.impl.length and value is not None:
            value = value[0:self.impl.length]
        return value


class Collector():

    def __init__(self, lastname=None, firstname=None, organization=None, email=None, collection_date=None):
        # Last name of the collector.
        self.lastname = lastname
        # First name of the collector.
        self.firstname = firstname
        self.organization = organization
        self.email = email
        self.collection_date = collection_date


class Coralvcf_allele():

    def __init__(self, chr=None, pos=None, sample_id=None, ref=None, alt=None, qual=None,
        filter=None, ac=None, an=None, bqb=None, dp=None, hob=None, icb=None, idf=None,
        imf=None, indel=None, mq=None, mqof=None, mqb=None, mqsb=None, rpb=None, sgb=None,
        vdb=None):
        self.chr = chr
        self.pos = pos
        self.sample_id = sample_id
        self.ref = ref
        self.alt = alt
        self.qual = qual
        self.filter = filter
        self.ac = ac
        self.an = an
        self.bqb = bqb
        self.dp = dp
        self.hob = hob
        self.icb = icb
        self.idf = idf
        self.imf = imf
        self.indel = indel
        self.mq = mq
        self.mqof = mqof
        self.mqb = mqb
        self.mqsb = mqsb
        self.rpb = rpb
        self.sgb = sgb
        self.vdb = vdb


class Experiment():

    def __init__(self, seq_facility=None, array_version=None, data_sharing=None, data_hold=None):
        # Description of experiment metadata.
        self.seq_facility = seq_facility
        self.array_version = array_version
        self.data_sharing = data_sharing
        self.data_hold = data_hold


class Genotype():

    def __init__(self, coral_mlg_clonal_id=None, symbiot_mlg_clonal_id=None, genetic_coral_species_call=None,
        percent_missing_data=None, percent_apalm=None, percent_acerv=None, percent_mixed=None):
        # Description of experiment metadata.
        self.coral_mlg_clonal_id = coral_mlg_clonal_id
        self.symbiot_mlg_clonal_id = symbiot_mlg_clonal_id
        self.genetic_coral_species_call = genetic_coral_species_call
        self.percent_missing_data = percent_missing_data
        self.percent_apalm = percent_apalm
        self.percent_acerv = percent_acerv
        self.percent_mixed = percent_mixed


class Person():

    def __init__(self, lastname=None, firstname=None, organization=None, email=None):
        # Last name of the collector.
        self.lastname = lastname
        # First name of the collector.
        self.firstname = firstname
        self.organization = organization
        self.email = email


class Reef():

    def __init__(self, name=None, region=None, latitude=None, longitude=None):
        self.name = name
        self.region = region
        # Latitude in decimal degrees, WGS84 datum, geographic location of reef.
        self.latitude = latitude
        # Longitude in decimal degrees, WGS84 datum, geographic location of reef.
        self.longitude = longitude


class Colony():

    def __init__(self, latitude=None, longitude=None, depth=None):
        # Latitude in decimal degrees, WGS84 datum, geographic location of
        # the wild donor colony from which the sample was collected.
        self.latitude = latitude
        # Longitude in decimal degrees, WGS84 datum, geographic location of
        # the wild donor colony from which the sample was collected.
        self.longitude = longitude
        # Depth in meters of the wild donor colony.
        self.depth = depth


class Sample():

    def __init__(self, sample_id=None, genotype_id=None, experiment_id=None, reef_id=None, colony_id=None,
        taxonomy_id=None, collector_id=None, collection_date=None, user_specimen_id=None, depth=None,
        dna_extraction_method=None, dna_concentration=None, duplicate_sample=None):
        self.sample_id = sample_id
        self.genotype_id = genotype_id
        self.experiment_id = experiment_id
        self.reef_id = reef_id
        self.colony_id = colony_id
        self.taxonomy_id = taxonomy_id
        self.collector_id = collector_id
        self.collection_date = collection_date
        # User specific identification for samples, should limit to some sort of length.
        self.user_specimin_id = user_specimin_id
        # Depth in meters of the sample for this entry.
        self.depth = depth
        # Method used to extract DNA.
        self.dna_extraction_method = dna_extraction_method
        # DNA concentration in ug.
        self.dna_concentration = dna_concentration
        # Indicate whether DNA from this colony or sample has been run before
        self.duplicate_sample = duplicate_sample


class Taxonomy():

    def __init__(self, species_name=None, genus_name=None):
        self.species_name = species_name
        self.genus_name = genus_name


Collector.table = Table("collector", metadata,
    Column("id", Integer, primary_key=True),
    Column("create_time", DateTime, default=now),
    Column("update_time", DateTime, default=now, onupdate=now),
    Column("person_id", Integer, ForeignKey("person.id"), index=True),
    Column("contact_id", Integer, ForeignKey("person.id"), index=True))


Colony.table = Table("colony", metadata,
    Column("id", Integer, primary_key=True),
    Column("latitude", Numeric(15, 10)),
    Column("longitude", Numeric(15, 10)),
    Column("depth", Integer))


Coralvcf_allele.table = Table("coralvcf_allele", metadata,
    Column("id", Integer, primary_key=True),
    Column("create_time", DateTime, default=now),
    Column("update_time", DateTime, default=now, onupdate=now),
    Column("chr", TrimmedString(255)),
    Column("pos", Integer),
    Column("sample_id", TrimmedString(255), index=True, nullable=False),
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


Experiment.table = Table("experiment", metadata,
    Column("id", Integer, primary_key=True),
    Column("create_time", DateTime, default=now),
    Column("update_time", DateTime, default=now, onupdate=now),
    Column("seq_facility", String),
    Column("array_version", TrimmedString(255)),
    Column("data_sharing", TrimmedString(255)),
    Column("data_hold", TrimmedString(255)))


Genotype.table = Table("genotype", metadata,
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


Person.table = Table("person", metadata,
    Column("id", Integer, primary_key=True),
    Column("create_time", DateTime, default=now),
    Column("update_time", DateTime, default=now, onupdate=now),
    Column("lastname", TrimmedString(255)),
    Column("firstname", TrimmedString(255)),
    Column("organization", TrimmedString(255)),
    Column("email", TrimmedString(255)))


Reef.table = Table("reef", metadata,
    Column("id", Integer, primary_key=True),
    Column("name", TrimmedString(255)),
    Column("region", TrimmedString(255)),
    Column("latitude", Numeric(15, 10)),
    Column("longitude", Numeric(15, 10)))


Sample.table = Table("sample", metadata,
    Column("id", Integer, primary_key=True),
    Column("create_time", DateTime, default=now),
    Column("update_time", DateTime, default=now, onupdate=now),
    Column("sample_id", TrimmedString(255), index=True, nullable=False),
    Column("genotype_id", Integer, ForeignKey("genotype.id"), index=True),
    Column("experiment_id", Integer, ForeignKey("experiment.id"), index=True),
    Column("reef_id", Integer, ForeignKey("reef.id"), index=True),
    Column("colony_id", Integer, ForeignKey("colony.id"), index=True),
    Column("taxonomy_id", Integer, ForeignKey("taxonomy.id"), index=True),
    Column("collector_id", Integer, ForeignKey("collector.id"), index=True),
    Column("collection_date", DateTime),
    Column("user_specimen_id", TrimmedString(255)),
    Column("depth", Integer),
    Column("dna_extraction_method", TrimmedString(255)),
    Column("dna_concentration", Numeric(10, 5)),
    Column("duplicate_sample", Boolean))


Taxonomy.table = Table("taxonomy", metadata,
    Column("id", Integer, primary_key=True),
    Column("create_time", DateTime, default=now),
    Column("update_time", DateTime, default=now, onupdate=now),
    Column("species_name", TrimmedString(255)),
    Column("genus_name", TrimmedString(255)))


# Now that the tables are defined, define the mappers
# and set up the relationships between the model objects.
#mapper(Sample, Sample.table, properties=dict(
#    colony=relation(Colony)
#)) 

metadata.create_all()
# conn = engine.connect()
# do more stuff here...


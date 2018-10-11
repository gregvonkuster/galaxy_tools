import datetime
import logging

from sqlalchemy import Boolean, Column, DateTime, ForeignKey, Integer, MetaData, Numeric, String, Table, TEXT

# Need our custom types, but don't import anything else from model
from galaxy.model.custom_types import JSONType, MetadataType, TrimmedString

now = datetime.datetime.utcnow
log = logging.getLogger(__name__)
metadata = MetaData()


Collector_table = Table("collector", metadata,
    Column("id", Integer, primary_key=True),
    Column("create_time", DateTime, default=now),
    Column("update_time", DateTime, default=now, onupdate=now),
    Column("person_id", Integer, ForeignKey("person.id"), index=True),
    Column("contact_id", Integer, ForeignKey("person.id"), index=True))


Colony_table = Table("colony", metadata,
    Column("id", Integer, primary_key=True),
    Column("create_time", DateTime, default=now),
    Column("update_time", DateTime, default=now, onupdate=now),
    Column("latitude", Numeric(15, 10)),
    Column("longitude", Numeric(15, 10)),
    Column("depth", Integer),
    Column("reef_id", Integer, ForeignKey("reef.id"), index=True))


Coralvcf_allele_table = Table("coralvcf_allele", metadata,
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


Experiment_table = Table("experiment", metadata,
    Column("id", Integer, primary_key=True),
    Column("create_time", DateTime, default=now),
    Column("update_time", DateTime, default=now, onupdate=now),
    Column("seq_facility", String),
    Column("array_version", TrimmedString(255)),
    Column("data_sharing", TrimmedString(255)),
    Column("data_hold", TrimmedString(255)))


Genotype_table = Table("genotype", metadata,
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


Idx_annotation_table = Table("idx_annotation", metadata,
    Column("id", Integer, primary_key=True),
    Column("create_time", DateTime, default=now),
    Column("update_time", DateTime, default=now, onupdate=now),
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


Person_table = Table("person", metadata,
    Column("id", Integer, primary_key=True),
    Column("create_time", DateTime, default=now),
    Column("update_time", DateTime, default=now, onupdate=now),
    Column("lastname", TrimmedString(255)),
    Column("firstname", TrimmedString(255)),
    Column("organization", TrimmedString(255)),
    Column("email", TrimmedString(255)))


Probe_annotation_table = Table("probe_annotation", metadata,
    Column("id", Integer, primary_key=True),
    Column("create_time", DateTime, default=now),
    Column("update_time", DateTime, default=now, onupdate=now),
    Column("probe_set_id", TrimmedString(255)),
    Column("affy_snp_id", TrimmedString(255)),
    Column("chr_id", Integer),
    Column("start", Integer),
    Column("stop", Integer),
    Column("strand", TrimmedString(255)),
    Column("dbsnp_rs_id", TrimmedString(255)),
    Column("dbsnp_loctype", Integer),
    Column("in_hapmap", TrimmedString(255)),
    Column("strand_vs_dbsnp", TrimmedString(255)),
    Column("probe_count", Integer),
    Column("cytoband", TrimmedString(255)),
    Column("chrx_par", Integer),
    Column("flank", TrimmedString(255)),
    Column("allele_a", TrimmedString(255)),
    Column("allele_b", TrimmedString(255)),
    Column("ref_allele", TrimmedString(255)),
    Column("alt_allele", TrimmedString(255)),
    Column("associated_gene", TrimmedString(255)),
    Column("genetic_map", TrimmedString(255)),
    Column("microsatellite", TrimmedString(255)),
    Column("heterozygous_allele_frequencies", TrimmedString(255)),
    Column("allele_frequency_count", TrimmedString(255)),
    Column("allele_frequencies", TrimmedString(255)),
    Column("minor_allele", TrimmedString(255)),
    Column("minor_allele_frequency", TrimmedString(255)),
    Column("omim", TrimmedString(255)),
    Column("biomedical", TrimmedString(255)),
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


Reef_table = Table("reef", metadata,
    Column("id", Integer, primary_key=True),
    Column("create_time", DateTime, default=now),
    Column("update_time", DateTime, default=now, onupdate=now),
    Column("name", TrimmedString(255)),
    Column("region", TrimmedString(255)),
    Column("latitude", Numeric(15, 10)),
    Column("longitude", Numeric(15, 10)))


Sample_table = Table("sample", metadata,
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
    Column("duplicate_sample", Boolean))


Taxonomy_table = Table("taxonomy", metadata,
    Column("id", Integer, primary_key=True),
    Column("create_time", DateTime, default=now),
    Column("update_time", DateTime, default=now, onupdate=now),
    Column("species_name", TrimmedString(255)),
    Column("genus_name", TrimmedString(255)))


def upgrade(migrate_engine):
    print(__doc__)
    metadata.bind = migrate_engine
    metadata.create_all()

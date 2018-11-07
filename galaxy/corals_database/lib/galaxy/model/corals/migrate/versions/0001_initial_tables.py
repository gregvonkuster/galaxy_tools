"""
Migration script to create the initial stag database and populate the
probe_annotation table from an external CSV file named probe_annotation.csv.
"""
from __future__ import print_function

import datetime
import logging
import os

from sqlalchemy import Boolean, Column, DateTime, ForeignKey, Integer, MetaData, Numeric, String, Table, TEXT

# Need our custom types, but don't import anything else from model
from galaxy.model.custom_types import JSONType, MetadataType, TrimmedString

now = datetime.datetime.utcnow
log = logging.getLogger(__name__)
metadata = MetaData()

# The current working direectory is the Galaxy
# installation root, so the following file must
# exist from that location.
GENERAL_SEED_DATA_FILE = "stag_database_seed_data/general_seed_data_file.csv"
PROBE_ANNOTATION_DATA_FILE = "stag_database_seed_data/probe_annotation.csv"

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
    Column("latitude", Numeric(15, 6)),
    Column("longitude", Numeric(15, 6)),
    Column("depth", Integer),
    Column("reef_id", Integer, ForeignKey("reef.id"), index=True))


Experiment_table = Table("experiment", metadata,
    Column("id", Integer, primary_key=True),
    Column("create_time", DateTime, default=now),
    Column("update_time", DateTime, default=now, onupdate=now),
    Column("seq_facility", String),
    Column("array_version", TrimmedString(255)),
    Column("data_sharing", TrimmedString(255)),
    Column("data_hold", TrimmedString(255)))


Fragment_table = Table("fragment", metadata,
    Column("id", Integer, primary_key=True),
    Column("create_time", DateTime, default=now),
    Column("update_time", DateTime, default=now, onupdate=now),
    Column("colony_id", Integer, ForeignKey("colony.id"), index=True))


Genotype_table = Table("genotype", metadata,
    Column("id", Integer, primary_key=True),
    Column("create_time", DateTime, default=now),
    Column("update_time", DateTime, default=now, onupdate=now),
    Column("coral_mlg_clonal_id", TrimmedString(255)),
    Column("symbiot_mlg_clonal_id", TrimmedString(255)),
    Column("genetic_coral_species_call", TrimmedString(255)),
    Column("percent_missing_data", Numeric(10, 6)),
    Column("percent_apalm", Numeric(10, 6)),
    Column("percent_acerv", Numeric(10, 6)),
    Column("percent_mixed", Numeric(10, 6)))


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


Reef_table = Table("reef", metadata,
    Column("id", Integer, primary_key=True),
    Column("create_time", DateTime, default=now),
    Column("update_time", DateTime, default=now, onupdate=now),
    Column("name", TrimmedString(255)),
    Column("region", TrimmedString(255)),
    Column("latitude", Numeric(15, 6)),
    Column("longitude", Numeric(15, 6)))


Sample_table = Table("sample", metadata,
    Column("id", Integer, primary_key=True),
    Column("create_time", DateTime, default=now),
    Column("update_time", DateTime, default=now, onupdate=now),
    Column("sample_id", TrimmedString(255), index=True, nullable=False),
    Column("genotype_id", Integer, ForeignKey("genotype.id"), index=True),
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
    Column("public", Boolean))


Taxonomy_table = Table("taxonomy", metadata,
    Column("id", Integer, primary_key=True),
    Column("create_time", DateTime, default=now),
    Column("update_time", DateTime, default=now, onupdate=now),
    Column("species_name", TrimmedString(255)),
    Column("genus_name", TrimmedString(255)))


def nextval(migrate_engine, table, col='id'):
    if migrate_engine.name in ['postgres', 'postgresql']:
        return "nextval('%s_%s_seq')" % (table, col)
    elif migrate_engine.name in ['mysql', 'sqlite']:
        return "null"
    else:
        raise Exception('Unable to convert data for unknown database type: %s' % migrate_engine.name)


def localtimestamp(migrate_engine):
    if migrate_engine.name in ['mysql', 'postgres', 'postgresql']:
        return "LOCALTIMESTAMP"
    elif migrate_engine.name == 'sqlite':
        return "current_date || ' ' || current_time"
    else:
        raise Exception('Unable to convert data for unknown database type: %s' % migrate_engine.name)


def get_latest_id(migrate_engine, table):
    result = migrate_engine.execute("select id from %s order by id desc" % table)
    row = result.fetchone()
    if row:
        return row[0]
    else:
        raise Exception('Unable to get the latest id in the %s table.' % table)

def load_probe_annotation_table(migrate_engine):
    # Columns:
    # probeset_id, affy_snp_id, chr_id, start, strand,
    # flank, allele_a, allele_b, allele_frequencies, annotation_notes, allele_count,
    # ordered_alleles, chrtype, custchr, custid, custpos,
    # organism, pconvert, recommendation, refstr, snppriority
    base_cmd = "INSERT INTO probe_annotation VALUES (%s, %s, %s, '%s', '%s', %s, %s, '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s')"

    with open(PROBE_ANNOTATION_DATA_FILE, "r") as fh:
        for i, line in enumerate(fh):
            if i == 0:
                # Skip the header.
                continue
            line = line.rstrip('\r\n')
            items = line.split(",")
            probeset_id = items[0]
            affy_snp_id = items[1]
            chr_id = items[2]
            start = items[3]
            strand = items[5]
            flank = items[13]
            allele_a = items[14]
            allele_b = items[15]
            allele_frequencies = items[23]
            annotation_notes = items[28]
            allele_count = items[29]
            ordered_alleles = items[30]
            chrtype = items[31]
            custchr = items[32]
            custid = items[33]
            custpos = items[34]
            organism = items[35]
            pconvert = items[36]
            recommendation = items[37]
            refstr = items[38]
            snppriority = items[39]
            cmd = base_cmd % (nextval(migrate_engine, 'probe_annotation'),
                              localtimestamp(migrate_engine),
                              localtimestamp(migrate_engine),
                              probeset_id,
                              affy_snp_id,
                              chr_id,
                              start,
                              strand,
                              flank,
                              allele_a,
                              allele_b,
                              allele_frequencies,
                              annotation_notes,
                              allele_count,
                              ordered_alleles,
                              chrtype,
                              custchr,
                              custid,
                              custpos,
                              organism,
                              pconvert,
                              recommendation,
                              refstr,
                              snppriority)
            migrate_engine.execute(cmd)
    print("Inserted %d rows into the probe_annotation table." % i)


def load_seed_data(migrate_engine):
    # Columns in general_seed_data_file.:
    # sample_id, date_entered_db, user_specimen_id, duplicate_sample, matching_samples,
    # field_call, bcoral_genet_id, bsym_genet_id, reef, region,
    # latitude, longitude, geographic_origin, sample_location, latitude_outplant,
    # longitude_outplant, depth, dist_shore, disease_resist, bleach_resist,
    # mortality, tle, spawning, collector, org,
    # collection_date, contact_email, seq_facility, array_version, data_sharing,
    # data_hold, coral_mlg_clonal_id, symbio_mlg_clonal_id, genetic_coral_species_call, percent_missing_data,
    # percent_apalm, percent_acerv, percent_mixed

    collector_table_inserts = 0
    experiment_table_inserts = 0
    person_table_inserts = 0
    reef_table_inserts = 0

    with open(GENERAL_SEED_DATA_FILE, "r") as fh:
        for i, line in enumerate(fh):
            if i == 0:
                # Skip the header.
                continue
            line = line.rstrip('\r\n')
            items = line.split(",")
            sample_id = items[0]
            data_entered_db = items[1]
            user_specimen_id = items[2]
            duplicate_sample = items[3]
            matching_samples = items[4]
            field_call = items[5]
            bcoral_genet_id = items[6]
            bsym_genet_id = items[7]
            reef = items[8]
            region = items[9]
            latitude = "%6f" % items[10]
            longitude = "%6f" % items[11]
            geographic_origin = items[12]
            sample_location = items[13]
            latitude_outplant = "%6f" % items[14]
            longitude_outplant = "%6f" % items[15]
            depth = items[16]
            dist_shore = items[17]
            disease_resist = items[18]
            bleach_resist = items[19]
            mortality = items[20]
            tle = items[21]
            spawning = items[22]
            collector = items[23]
            org = items[24]
            collection_date = items[25]
            contact_email = items[26].lower()
            seq_facility = items[27]
            array_version = items[28]
            data_sharing = items[29]
            data_hold = items[30]
            coral_mlg_clonal_id = items[31]
            symbio_mlg_clonal_id = items[32]
            genetic_coral_species_call = items[33]
            percent_missing_data = "%6f" % items[34]
            percent_apalm = "%6f" % items[35]
            percent_acerv = "%6f" % items[36]
            percent_mixed = "%6f" % items[37]
            # See if we need to add a row to the person table.
            if collector.find(" ") > 0:
                # We have a first and last name spearated by a space.
                first_last = collector.split(" ")
                first_name = first_last[0]
                last_name = first_last[1]
                cmd = "SELECT id FROM person WHERE last_name = '%s' AND first_name = '%s' AND email = '%s'" % (last_name, first_name, email)
            else:
                # We have a last name with no first name.
                cmd = "SELECT id FROM person WHERE last_name = '%s' and email = '%s'" % (last_name, email)
                result = migrate_engine.execute(cmd)
                row = result.fetchone()
                if row:
                    person_id = row[0]
                else:
                    # Add a row to the person table.
                    cmd = "INSERT INTO person VALUES (%s, %s, %s, '%s', '%s', '%s', '%s')"
                    cmd = cmd % (nextval(migrate_engine, 'person'),
                                 localtimestamp(migrate_engine),
                                 localtimestamp(migrate_engine),
                                 last_name,
                                 first_name,
                                 org,
                                 contact_email)
                    migrate_engine.execute(cmd)
                    person_table_inserts += 1
                    person_id = get_latest_id(migrate_engine, "person")
            # Add a row to the collector table.
            cmd = "INSERT INTO collector VALUES (%s, %s, %s, %s, %s)"
            cmd = cmd % (nextval(migrate_engine, 'person'),
                         localtimestamp(migrate_engine),
                         localtimestamp(migrate_engine),
                         person_id,
                         pserson_id)
            migrate_engine.execute(cmd)
            collector_table_inserts += 1
            # See if we need to add a row to the experiment table.
Experiment_table = Table("experiment", metadata,
    Column("id", Integer, primary_key=True),
    Column("create_time", DateTime, default=now),
    Column("update_time", DateTime, default=now, onupdate=now),
    Column("seq_facility", String),
    Column("array_version", TrimmedString(255)),
    Column("data_sharing", TrimmedString(255)),
    Column("data_hold", TrimmedString(255)))
            cmd = "SELECT id FROM experiment WHERE seq_facility = '%s' AND array_version = '%s' AND data_sharing = '%s' AND data_hold = '%s'"
            cmd = cmd % (seq_facility, array_version, data_sharing, data_hold)
            result = migrate_engine.execute(cmd)
            row = result.fetchone()
            if row:
                experiment_id = row[0]
            else:
                # Add a row to the experiment table.
                cmd = "INSERT INTO experiment VALUES (%s, %s, %s, '%s', '%s', '%s', '%s')"
                cmd = cmd % (nextval(migrate_engine, 'experiment'),
                             localtimestamp(migrate_engine),
                             localtimestamp(migrate_engine),
                             seq_facility,
                             array_version,
                             data_sharing,
                             data_hold)
                migrate_engine.execute(cmd)
                experiment_table_inserts += 1
                experiment_id = get_latest_id(migrate_engine, "experiment")
            # See if we need to add a row to the reef table.
            cmd = "SELECT id FROM reef WHERE name = '%s'" % reef
            result = migrate_engine.execute(cmd)
            row = result.fetchone()
            if row:
                reef_id = row[0]
            else:
                # Add a row to the person table.
                cmd = "INSERT INTO reef VALUES (%s, %s, %s, '%s', '%s', %s, %s)"
                cmd = cmd % (nextval(migrate_engine, 'reef'),
                             localtimestamp(migrate_engine),
                             localtimestamp(migrate_engine),
                             reef,
                             region,
                             latitude,
                             longitude)
                migrate_engine.execute(cmd)
                reef_table_inserts += 1
                reef_id = get_latest_id(migrate_engine, "reef")

    print("Inserted %d rows into the collector table." % collector_table_inserts)
    print("Inserted %d rows into the experiment table." % experiment_table_inserts)
    print("Inserted %d rows into the person table." % person_table_inserts)
    print("Inserted %d rows into the reef table." % reef_table_inserts)
            
def downgrade(migrate_engine):
    pass

def upgrade(migrate_engine):
    print(__doc__)
    metadata.bind = migrate_engine
    metadata.create_all()
    load_seed_data(migrate_engine)
    load_probe_annotation_table(migrate_engine)

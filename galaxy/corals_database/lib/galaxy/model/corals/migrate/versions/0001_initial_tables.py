"""
Migration script to create the initial stag database and populate the database
with production seed data from a file named general_seed_data_file.tabular.  The
probe_annotation table is also populated from an external CSV file named
probe_annotation.csv.
"""
from __future__ import print_function

import datetime
import dateutil.parser
import logging
import time

from sqlalchemy import Boolean, Column, DateTime, ForeignKey, Integer, MetaData, Numeric, String, Table, TIMESTAMP
from sqlalchemy import sql

# Need our custom types, but don't import anything else from model
from galaxy.model.custom_types import TrimmedString

now = datetime.datetime.utcnow
log = logging.getLogger(__name__)
metadata = MetaData()

# Get current date plus one year for possible insertion
# into the public_after_date column of the sample table.
# The default behavior is for the value of the public
# column to be True and the public_after_date to be NULL,
# making the sample "public".  HOwever, the user can
# set the value of the public column to False and optionally
# set a date after which the sample becomes public in the
# Affymetrix 96 well plate metadata file associated with
# the sample.  If the value of the public column is set
# to False, but no date is set, the default date will be 1
# year from the time the row is inserted into the table. 
today = datetime.date.today()
try:
    # Return the same day of the year.
    year_from_now = today.replace(year=today.year+1)
except Exception:
    # Handle leap years.
    year_from_now = today + (datetime.date(today.year + 1, 1, 1) - datetime.date(today.year, 1, 1))

# The current working direectory is the Galaxy
# installation root, so the following file must
# exist from that location.
GENERAL_SEED_DATA_FILE = "stag_database_seed_data/general_seed_data_file.tabular"
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
    Column("array_version", TrimmedString(255)))


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
    Column("coral_mlg_rep_sample_id", TrimmedString(255)),
    Column("symbio_mlg_clonal_id", TrimmedString(255)),
    Column("symbio_mlg_rep_sample_id", TrimmedString(255)),
    Column("genetic_coral_species_call", TrimmedString(255)))


Person_table = Table("person", metadata,
    Column("id", Integer, primary_key=True),
    Column("create_time", DateTime, default=now),
    Column("update_time", DateTime, default=now, onupdate=now),
    Column("last_name", TrimmedString(255)),
    Column("first_name", TrimmedString(255)),
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


Phenotype_table = Table("phenotype", metadata,
    Column("id", Integer, primary_key=True),
    Column("create_time", DateTime, default=now),
    Column("update_time", DateTime, default=now, onupdate=now),
    Column("disease_resist", TrimmedString(255)),
    Column("bleach_resist", TrimmedString(255)),
    Column("mortality", TrimmedString(255)),
    Column("tle", TrimmedString(255)),
    Column("spawning", Boolean),
    Column("sperm_motility", Numeric(15, 6), nullable=False),
    Column("healing_time", Numeric(15, 6), nullable=False))


Sample_table = Table("sample", metadata,
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
    Column("percent_heterozygous_coral", Numeric(15, 6)),
    Column("percent_heterozygous_sym", Numeric(15, 6)))


Taxonomy_table = Table("taxonomy", metadata,
    Column("id", Integer, primary_key=True),
    Column("create_time", DateTime, default=now),
    Column("update_time", DateTime, default=now, onupdate=now),
    Column("species_name", TrimmedString(255)),
    Column("genus_name", TrimmedString(255)))


def boolean(migrate_engine, value):
    if migrate_engine.name in ['mysql', 'postgres', 'postgresql']:
        return value
    elif migrate_engine.name == 'sqlite':
        if value in ['True', 'true']:
            return 1
        return 0
    else:
        raise Exception('Unable to convert data for unknown database type: %s' % migrate_engine.name)


def string_as_bool(string):
    if str(string).lower() in ('true', 'yes', 'on', '1'):
        return True
    else:
        return False


def convert_date_string_for_database(date_string):
    # The value of date_string is %y/%m/%d with
    # the year being 2 digits (yikes!).
    fixed_century = "20%s" % date_string
    fixed_date = fixed_century.replace("/", "-")
    # Convert the string to a format required for
    # inserting into the database.
    database_format = dateutil.parser.parse(fixed_date)
    return str(database_format)


def get_latest_id(migrate_engine, table):
    result = migrate_engine.execute("select id from %s order by id desc" % table)
    row = result.fetchone()
    if row:
        return row[0]
    else:
        raise Exception('Unable to get the latest id in the %s table.' % table)


def get_primary_id(migrate_engine, table, cmd):
    result = migrate_engine.execute(cmd)
    row = result.fetchone()
    if row:
        return row[0]
    else:
        return None


def localtimestamp(migrate_engine):
    if migrate_engine.name in ['mysql', 'postgres', 'postgresql']:
        return "LOCALTIMESTAMP"
    elif migrate_engine.name == 'sqlite':
        return "current_date || ' ' || current_time"
    else:
        raise Exception('Unable to convert data for unknown database type: %s' % migrate_engine.name)


def nextval(migrate_engine, table, col='id'):
    if migrate_engine.name in ['postgres', 'postgresql']:
        return "nextval('%s_%s_seq')" % (table, col)
    elif migrate_engine.name in ['mysql', 'sqlite']:
        return "null"
    else:
        raise Exception('Unable to convert data for unknown database type: %s' % migrate_engine.name)


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
    # user_specimen_id, field_call, bcoral_genet_id, bsym_genet_id, reef,
    # region, latitude, longitude, geographic_origin, colony_location,
    # latitude_outplant, longitude_outplant, depth, dist_shore, disease_resist,
    # bleach_resist, mortality, tle, spawning, collector, org,
    # collection_date, contact_email, seq_facility, array_version, public,
    # public_after_date, coral_mlg_clonal_id, symbio_mlg_clonal_id, genetic_coral_species_call, percent_missing_data,
    # percent_apalm, percent_acerv, percent_mixed, affy_id

    collector_table_inserts = 0
    colony_table_inserts = 0
    experiment_table_inserts = 0
    genotype_table_inserts = 0
    person_table_inserts = 0
    phenotype_table_inserts = 0
    reef_table_inserts = 0
    sample_table_inserts = 0
    SAMPLE_ID = 10000
    # The following are not included in the seed data,
    # so negative values are used since they are invalid
    # and so can act as so-called NULL values.
    SPERM_MOTILITY = -9.0
    HEALING_TIME = -9.0

    with open(GENERAL_SEED_DATA_FILE, "r") as fh:
        for i, line in enumerate(fh):
            if i == 0:
                # Skip the header.
                continue
            line = line.rstrip('\r\n')
            items = line.split("\t")
            # Automatically generate the sample_id.
            sample_id = "A%d" % SAMPLE_ID
            SAMPLE_ID += 1
            user_specimen_id = items[0]
            if len(items[1]) == 0:
                field_call = sql.null()
            else:
                field_call = items[1]
            if len(items[2]) == 0:
                bcoral_genet_id = sql.null()
            else:
                bcoral_genet_id = items[2]
            if len(items[3]) == 0:
                bsym_genet_id = sql.null()
            else:
                bsym_genet_id = items[3]
            reef = items[4]
            region = items[5]
            try:
                latitude = "%6f" % float(items[6])
            except Exception:
                latitude = sql.null()
            try:
                longitude = "%6f" % float(items[7])
            except Exception:
                longitude = sql.null()
            if len(items[8]) == 0:
                geographic_origin = sql.null()
            else:
                geographic_origin = items[8]
            if len(items[9]) == 0:
                colony_location = sql.null()
            else:
                colony_location = items[9]
            try:
                latitude_outplant = "%6f" % float(items[10])
            except Exception:
                latitude_outplant = sql.null()
            try:
                longitude_outplant = "%6f" % float(items[11])
            except Exception:
                longitude_outplant = sql.null()
            try:
                depth = int(items[12])
            except Exception:
                depth = 0
            if len(items[13]) == 0:
                dist_shore = sql.null()
            else:
                dist_shore = items[13]
            disease_resist = items[14]
            bleach_resist = items[15]
            mortality = items[16]
            tle = items[17]
            # Convert original spawning value to Boolean.
            spawning = string_as_bool(items[18])
            collector = items[19]
            org = items[20]
            try:
                collection_date = convert_date_string_for_database(items[21])
            except Exception:
                collection_date = localtimestamp(migrate_engine)
            contact_email = items[22]
            seq_facility = items[23]
            array_version = items[24]
            # Convert original public value to Boolean.
            public = string_as_bool(items[25])
            if public:
                public_after_date = sql.null()
            else:
                if len(items[26]) == 0:
                    # Set the value of public_after_date to the default.
                    public_after_date = year_from_now
                else:
                    public_after_date = convert_date_string_for_database(items[26])
            coral_mlg_clonal_id = items[27]
            symbio_mlg_clonal_id = items[28]
            genetic_coral_species_call = items[29]
            try:
                percent_missing_data = "%6f" % float(items[30])
            except Exception:
                percent_missing_data = sql.null()
            try:
                percent_apalm = "%6f" % float(items[31])
            except Exception:
                percent_apalm = sql.null()
            try:
                percent_acerv = "%6f" % float(items[32])
            except Exception:
                percent_acerv = sql.null()
            try:
                percent_mixed = "%6f" % float(items[33])
            except Exception:
                percent_mixed = sql.null()
            affy_id = items[34]

            # Process the experiment items.  Dependent tables: sample.
            table = "experiment"
            # See if we need to add a row to the experiment table.
            cmd = "SELECT id FROM experiment WHERE seq_facility = '%s' AND array_version = '%s'"
            cmd = cmd % (seq_facility, array_version)
            experiment_id = get_primary_id(migrate_engine, table, cmd)
            if experiment_id is None:
                # Add a row to the experiment table.
                cmd = "INSERT INTO experiment VALUES (%s, %s, %s, '%s', '%s')"
                cmd = cmd % (nextval(migrate_engine, table),
                             localtimestamp(migrate_engine),
                             localtimestamp(migrate_engine),
                             seq_facility,
                             array_version)
                migrate_engine.execute(cmd)
                experiment_table_inserts += 1
                experiment_id = get_latest_id(migrate_engine, table)

            # Process the genotype items.  Dependent tables: sample.
            table = "genotype"
            # See if we need to add a row to the table.
            cmd = "SELECT id FROM genotype WHERE coral_mlg_clonal_id = '%s' AND symbio_mlg_clonal_id = '%s' AND genetic_coral_species_call = '%s'"
            cmd = cmd % (coral_mlg_clonal_id,
                         symbio_mlg_clonal_id,
                         genetic_coral_species_call)
            genotype_id = get_primary_id(migrate_engine, table, cmd)
            if genotype_id is None:
                # Add a row to the table.
                cmd = "INSERT INTO genotype VALUES (%s, %s, %s, '%s', '%s', '%s')"
                cmd = cmd % (nextval(migrate_engine, table),
                             localtimestamp(migrate_engine),
                             localtimestamp(migrate_engine),
                             coral_mlg_clonal_id,
                             symbio_mlg_clonal_id,
                             genetic_coral_species_call)
                migrate_engine.execute(cmd)
                genotype_table_inserts += 1
                genotype_id = get_latest_id(migrate_engine, table)

            # Process the phenotype items.  Dependent tables: sample.
            table = "phenotype"
            # See if we need to add a row to the table.
            cmd = "SELECT id FROM phenotype WHERE disease_resist = '%s' AND bleach_resist = '%s' AND mortality = '%s'"
            cmd += " AND tle = '%s' AND spawning = '%s'"
            cmd = cmd % (disease_resist,
                         bleach_resist,
                         mortality,
                         tle,
                         spawning)
            phenotype_id = get_primary_id(migrate_engine, table, cmd)
            if phenotype_id is None:
                # Add a row to the table.
                cmd = "INSERT INTO phenotype VALUES (%s, %s, %s, '%s', '%s', '%s', '%s', '%s', %s, %s)"
                cmd = cmd % (nextval(migrate_engine, table),
                             localtimestamp(migrate_engine),
                             localtimestamp(migrate_engine),
                             disease_resist,
                             bleach_resist,
                             mortality,
                             tle,
                             spawning,
                             SPERM_MOTILITY,
                             HEALING_TIME)
                migrate_engine.execute(cmd)
                phenotype_table_inserts += 1
                phenotype_id = get_latest_id(migrate_engine, table)

            # Process the person items.  Dependent tables: collector.
            table = "person"
            # See if we need to add a row to the table.
            if collector.find(" ") > 0:
                # We have a first and last name spearated by a space.
                first_last = collector.split(" ")
                first_name = first_last[0]
                last_name = first_last[1]
                cmd = "SELECT id FROM person WHERE last_name = '%s' AND first_name = '%s' AND email = '%s'" % (last_name, first_name, contact_email)
            else:
                # We have a last name with no first name.
                if len(collector) > 0:
                    last_name = collector
                else:
                    last_name = 'Unknown'
                first_name = sql.null()
                cmd = "SELECT id FROM person WHERE last_name = '%s' and email = '%s'" % (last_name, contact_email)
            person_id = get_primary_id(migrate_engine, table, cmd)
            if person_id is None:
                # Add a row to the table.
                cmd = "INSERT INTO person VALUES (%s, %s, %s, '%s', '%s', '%s', '%s')"
                cmd = cmd % (nextval(migrate_engine, table),
                             localtimestamp(migrate_engine),
                             localtimestamp(migrate_engine),
                             last_name,
                             first_name,
                             org,
                             contact_email)
                migrate_engine.execute(cmd)
                person_table_inserts += 1
                person_id = get_latest_id(migrate_engine, table)

            # Process the collector items.  Dependent tables: sample.
            table = "collector"
            # See if we need to add a row to the table.
            cmd = "SELECT id FROM collector WHERE person_id = %s" % person_id
            collector_id = get_primary_id(migrate_engine, table, cmd)
            if collector_id is None:
                # Add a row to the table.
                cmd = "INSERT INTO collector VALUES (%s, %s, %s, %s, %s)"
                cmd = cmd % (nextval(migrate_engine, table),
                             localtimestamp(migrate_engine),
                             localtimestamp(migrate_engine),
                             person_id,
                             person_id)
                migrate_engine.execute(cmd)
                collector_table_inserts += 1
                collector_id = get_latest_id(migrate_engine, table)

            # Process the reef items.  Dependent tables: colony.
            table = "reef"
            # See if we need to add a row to the table.
            cmd = "SELECT id FROM reef WHERE name = '%s'" % reef
            reef_id = get_primary_id(migrate_engine, table, cmd)
            if reef_id is None:
                # Add a row to the table.
                cmd = "INSERT INTO reef VALUES (%s, %s, %s, '%s', '%s', %s, %s)"
                cmd = cmd % (nextval(migrate_engine, table),
                             localtimestamp(migrate_engine),
                             localtimestamp(migrate_engine),
                             reef,
                             region,
                             latitude,
                             longitude)
                migrate_engine.execute(cmd)
                reef_table_inserts += 1
                reef_id = get_latest_id(migrate_engine, table)

            # Process the colony items.  Dependent tables: fragment, sample.
            table = "colony"
            # See if we need to add a row to the table.
            cmd = "SELECT id FROM colony WHERE latitude = %s AND longitude = %s and reef_id = %s"
            cmd = cmd % (latitude_outplant,
                         longitude_outplant,
                         reef_id)
            colony_id = get_primary_id(migrate_engine, table, cmd)
            if colony_id is None:
                # Add a row to the table.
                cmd = "INSERT INTO colony VALUES (%s, %s, %s, %s, %s, %s, %s)"
                cmd = cmd % (nextval(migrate_engine, table),
                             localtimestamp(migrate_engine),
                             localtimestamp(migrate_engine),
                             latitude_outplant,
                             longitude_outplant,
                             depth,
                             reef_id)
                migrate_engine.execute(cmd)
                colony_table_inserts += 1
                colony_id = get_latest_id(migrate_engine, table)

            # Process the sample items.  Dependent tables: None.
            table = "sample"
            # See if we need to add a row to the table.
            cmd = "SELECT id FROM sample WHERE sample_id = '%s'" % sample_id
            sample_id_db = get_primary_id(migrate_engine, table, cmd)
            if sample_id_db is None:
                # Add a row to the table.  Values for
                # the following are not in the seed data.
                fragment_id = sql.null()
                taxonomy_id = sql.null()
                dna_extraction_method = sql.null()
                dna_concentration = sql.null()
                percent_missing_data_coral = sql.null()
                percent_missing_data_sym = sql.null()
                percent_reference_coral = sql.null()
                percent_reference_sym = sql.null()
                percent_alternative_coral = sql.null()
                percent_alternative_sym = sql.null()
                percent_heterozygous_coral = sql.null()
                percent_heterozygous_sym = sql.null()
                # id, create_time, update_time, affy_id, sample_id,
                # genotype_id, phenotype_id, experiment_id, colony_id, colony_location,
                # fragment_id, taxonomy_id, collector_id
                cmd = "INSERT INTO sample VALUES (%s, %s, %s, '%s', '%s', %s, %s, %s, %s, %s, %s, %s, %s, "
                if collection_date == "LOCALTIMESTAMP":
                    # collection_date
                    cmd += "%s, "
                else:
                    # collection_date
                    cmd += "'%s', "
                # user_specimen_id, depth, dna_extraction_method, dna_concentration, public,
                # public_after_date, percent_missing_data_coral, percent_missing_data_sym,
                # percent_reference_coral, percent_reference_sym, percent_alternative_coral,
                # percent_alternative_sym, percent_heterozygous_coral, percent_heterozygous_sym
                cmd += "'%s', %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)"
                cmd = cmd % (nextval(migrate_engine, table),
                             localtimestamp(migrate_engine),
                             localtimestamp(migrate_engine),
                             affy_id,
                             sample_id,
                             genotype_id,
                             phenotype_id,
                             experiment_id,
                             colony_id,
                             colony_location,
                             fragment_id,
                             taxonomy_id,
                             collector_id,
                             collection_date,
                             user_specimen_id,
                             depth,
                             dna_extraction_method,
                             dna_concentration,
                             public,
                             public_after_date,
                             percent_missing_data_coral,
                             percent_missing_data_sym,
                             percent_reference_coral,
                             percent_reference_sym,
                             percent_alternative_coral,
                             percent_alternative_sym,
                             percent_heterozygous_coral,
                             percent_heterozygous_sym)
                migrate_engine.execute(cmd)
                sample_table_inserts += 1
                sample_id = get_latest_id(migrate_engine, table)

    print("Inserted %d rows into the collector table." % collector_table_inserts)
    print("Inserted %d rows into the colony table." % colony_table_inserts)
    print("Inserted %d rows into the experiment table." % experiment_table_inserts)
    print("Inserted %d rows into the genotype table." % genotype_table_inserts)
    print("Inserted %d rows into the person table." % person_table_inserts)
    print("Inserted %d rows into the phenotype table." % phenotype_table_inserts)
    print("Inserted %d rows into the reef table." % reef_table_inserts)
    print("Inserted %d rows into the sample table." % sample_table_inserts)


def downgrade(migrate_engine):
    pass


def upgrade(migrate_engine):
    print(__doc__)
    metadata.bind = migrate_engine
    metadata.create_all()
    load_seed_data(migrate_engine)
    load_probe_annotation_table(migrate_engine)

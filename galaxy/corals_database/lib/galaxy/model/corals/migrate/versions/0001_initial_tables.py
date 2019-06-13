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

from sqlalchemy import (
    Boolean,
    Column,
    DateTime,
    ForeignKey,
    Integer,
    MetaData,
    Numeric,
    String,
    Table,
    Text
)
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
    year = today.year + 1
    year_from_now = today.replace(year=year)
except Exception:
    # Handle leap years.
    year_from_now = today + (datetime.date(today.year + 1, 1, 1) - datetime.date(today.year, 1, 1))

# The current working directory is the Galaxy
# installation root, so the following file must
# exist from that location.
ALLELES_SEED_DATA_FILE = "stag_database_seed_data/alleles_seed_data_file.tabular"
GENERAL_SEED_DATA_FILE = "stag_database_seed_data/general_seed_data_file.tabular"
PROBE_ANNOTATION_DATA_FILE = "stag_database_seed_data/probe_annotation.csv"
DEFAULT_MISSING_NUMERIC_VALUE = -9.000000

Allele_table = Table("allele", metadata,
                     Column("id", Integer, primary_key=True),
                     Column("create_time", DateTime, default=now),
                     Column("update_time", DateTime, default=now, onupdate=now),
                     Column("allele", Text))


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
                         Column("result_folder_name", TrimmedString(255)),
                         Column("plate_barcode", TrimmedString(255)))


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
                       Column("genetic_coral_species_call", TrimmedString(255)),
                       Column("bcoral_genet_id", TrimmedString(255)))


Person_table = Table("person", metadata,
                     Column("id", Integer, primary_key=True),
                     Column("create_time", DateTime, default=now),
                     Column("update_time", DateTime, default=now, onupdate=now),
                     Column("last_name", TrimmedString(255)),
                     Column("first_name", TrimmedString(255)),
                     Column("organization", TrimmedString(255)),
                     Column("email", TrimmedString(255)))


Phenotype_table = Table("phenotype", metadata,
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
                   Column("longitude", Numeric(15, 6)),
                   Column("geographic_origin", TrimmedString(255)))


Sample_table = Table("sample", metadata,
                     Column("id", Integer, primary_key=True),
                     Column("create_time", DateTime, default=now),
                     Column("update_time", DateTime, default=now, onupdate=now),
                     Column("affy_id", TrimmedString(255), index=True, nullable=False),
                     Column("sample_id", TrimmedString(255), index=True, nullable=False),
                     Column("allele_id", Integer, ForeignKey("allele.id"), index=True),
                     Column("genotype_id", Integer, ForeignKey("genotype.id"), index=True),
                     Column("phenotype_id", Integer, ForeignKey("phenotype.id"), index=True),
                     Column("experiment_id", Integer, ForeignKey("experiment.id"), index=True),
                     Column("colony_id", Integer, ForeignKey("colony.id"), index=True),
                     Column("colony_location", TrimmedString(255)),
                     Column("taxonomy_id", Integer, ForeignKey("taxonomy.id"), index=True),
                     Column("collector_id", Integer, ForeignKey("person.id"), index=True),
                     Column("collection_date", DateTime),
                     Column("user_specimen_id", TrimmedString(255)),
                     Column("registry_id", TrimmedString(255)),
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
                     Column("percent_heterozygous_sym", Numeric(15, 6)),
                     Column("field_call", TrimmedString(255)))


Symbio_genotype_table = Table("symbio_genotype", metadata,
                             Column("id", Integer, primary_key=True),
                             Column("create_time", DateTime, default=now),
                             Column("update_time", DateTime, default=now, onupdate=now),
                             Column("symbio_mlg_clonal_id", TrimmedString(255)),
                             Column("symbio_mlg_rep_sample_id", TrimmedString(255)),
                             Column("bsym_genet_id", TrimmedString(255)))


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


def handle_column_value(val, default=''):
    # Regarding the default value, a NULL value indicates n unknown value
    # and typically should not be confused with an empty string. Our application
    # does not need the concept of unknown value, so most columns are
    # non-nullable and our default is an empty string.
    if val in ["", "NA", "NULL"]:
        return default
    return val


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


def load_alleles_seed_data(migrate_engine):
    # Columns in alleles_seed_data_file.:
    # [0]user_specimen_id [1]alleles

    allele_table_inserts = 0
    sample_table_updates = 0

    with open(ALLELES_SEED_DATA_FILE, "r") as fh:
        for i, line in enumerate(fh):
            if i == 0:
                # Skip the header.
                continue
            line = line.rstrip('\r\n')
            items = line.split("\t")
            user_specimen_id = items[0]
            alleles = items[1]

            # Process the allele items.  Dependent tables: sample.
            table = "allele"
            # Add a row to the allele table.
            cmd = "INSERT INTO allele VALUES (%s, %s, %s, '%s')"
            cmd = cmd % (nextval(migrate_engine, table),
                         localtimestamp(migrate_engine),
                         localtimestamp(migrate_engine),
                         alleles)
            migrate_engine.execute(cmd)
            allele_table_inserts += 1
            allele_id = get_latest_id(migrate_engine, table)

            # Update the row in the sample table that contains
            # the user_specimen_id with the allele_id.
            cmd = "UPDATE sample SET allele_id = %s where user_specimen_id = '%s'" % (allele_id, user_specimen_id)
            migrate_engine.execute(cmd)
            sample_table_updates += 1

    print("Inserted %d rows into the allele table." % allele_table_inserts)
    print("Updated %d rows in the sample table." % sample_table_updates)


def load_general_seed_data(migrate_engine):
    # Columns in general_seed_data_file.:
    # [0]user_specimen_id [1]field_call [2]bcoral_genet_id [3]bsym_genet_id [4]reef
    # [5]region [6]latitude [7]longitude [8]geographic_origin [9]colony_location
    # [10]latitude_outplant [11]longitude_outplant [12]depth [13]disease_resist [14]bleach_resist
    # [15]mortality [16]tle [17]spawning [18]collector_last_name [19]collector_first_name
    # [20]organization [21]collection_date [22]email [23]seq_facility [24]array_version
    # [25]public [26]public_after_date [27]coral_mlg_clonal_id [28]symbio_mlg_clonal_id [29]genetic_coral_species_call
    # [30]percent_missing_data_coral [31]percent_reference_coral [32]percent_alternative_sym [33]percent_heterozygous_coral [34]affy_id
    # [35]coral_mlg_rep_sample_id [36]genus_name [37]species [38]sperm_motility [39]healing_time
    # [40]dna_extraction_method [41]dna_concentration [42]registry_id

    colony_table_inserts = 0
    experiment_table_inserts = 0
    genotype_table_inserts = 0
    person_table_inserts = 0
    phenotype_table_inserts = 0
    reef_table_inserts = 0
    sample_table_inserts = 0
    taxonomy_table_inserts = 0
    SAMPLE_ID = 10000

    # The following values are not included in the seed data.
    result_folder_name = ''
    plate_barcode = ''

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
                field_call = ''
            else:
                field_call = handle_column_value(items[1])
            if len(items[2]) == 0:
                bcoral_genet_id = ''
            else:
                bcoral_genet_id = handle_column_value(items[2])
            if len(items[3]) == 0:
                bsym_genet_id = ''
            else:
                bsym_genet_id = handle_column_value(items[3])
            reef = handle_column_value(items[4])
            region = handle_column_value(items[5])
            try:
                latitude = "%6f" % float(items[6])
            except Exception:
                # Latitude of Mueller Lab building.
                latitude = 40.799064
            try:
                # Longitude of Mueller Lab building.
                longitude = "%6f" % float(items[7])
            except Exception:
                longitude = -77.864431
            if len(items[8]) == 0:
                # TODO: find out if the default value for
                # geographihc origin should be "reef".
                geographic_origin = 'reef'
            else:
                geographic_origin = handle_column_value(items[8], default='reef')
                geographic_origin = geographic_origin.lower()
            if len(items[9]) == 0:
                colony_location = ''
            else:
                colony_location = items[9]
            try:
                depth = int(items[12])
            except Exception:
                depth = 0
            disease_resist = items[13]
            bleach_resist = items[14]
            mortality = items[15]
            tle = items[16]
            # Convert original spawning value to Boolean.
            spawning = string_as_bool(items[17])
            collector_last_name = handle_column_value(items[18])
            collector_first_name = handle_column_value(items[19])
            organization = handle_column_value(items[20])
            try:
                collection_date = convert_date_string_for_database(items[21])
            except Exception:
                collection_date = localtimestamp(migrate_engine)
            email = handle_column_value(items[22])
            seq_facility = items[23]
            array_version = handle_column_value(items[24])
            # Convert original public value to Boolean.
            public = string_as_bool(items[25])
            if public:
                public_after_date = 'NOW()'
            else:
                if len(items[26]) == 0:
                    # Set the value of public_after_date to the default.
                    public_after_date = year_from_now
                else:
                    public_after_date = convert_date_string_for_database(items[26])
            coral_mlg_clonal_id = handle_column_value(items[27])
            symbio_mlg_clonal_id = handle_column_value(items[28])
            genetic_coral_species_call = handle_column_value(items[29])
            try:
                percent_missing_data_coral = "%6f" % float(items[30])
            except Exception:
                percent_missing_data_coral = DEFAULT_MISSING_NUMERIC_VALUE
            try:
                percent_reference_coral = "%6f" % float(items[31])
            except Exception:
                percent_reference_coral = DEFAULT_MISSING_NUMERIC_VALUE
            try:
                percent_alternative_sym = "%6f" % float(items[32])
            except Exception:
                percent_alternative_sym = DEFAULT_MISSING_NUMERIC_VALUE
            try:
                percent_heterozygous_coral = "%6f" % float(items[33])
            except Exception:
                percent_heterozygous_coral = DEFAULT_MISSING_NUMERIC_VALUE
            affy_id = items[34]
            if len(items[35]) == 0:
                coral_mlg_rep_sample_id = ''
            else:
                coral_mlg_rep_sample_id = handle_column_value(items[35])
            if len(items[36]) == 0:
                genus_name = ''
            else:
                genus_name = handle_column_value(items[36])
            if len(items[37]) == 0:
                species_name = ''
            else:
                species_name = handle_column_value(items[37])
            try:
                sperm_motility = "%6f" % float(items[38])
            except Exception:
                sperm_motility = ''
            try:
                healing_time = "%6f" % float(items[39])
            except Exception:
                healing_time = ''
            if len(items[40]) == 0:
                dna_extraction_method = ''
            else:
                dna_extraction_method = handle_column_value(items[40])
            try:
                dna_concentration = "%6f" % float(items[41])
            except Exception:
                dna_concentration = ''
            if len(items[42]) == 0:
                registry_id = ''
            else:
                registry_id = handle_column_value(items[42])

            # Process the taxonomy items.  Dependent tables: sample.
            table = "taxonomy"
            # See if we need to add a row to the taxonomy table.
            cmd = "SELECT id FROM taxonomy WHERE species_name = '%s' AND genus_name = '%s'"
            cmd = cmd % (species_name, genus_name)
            taxonomy_id = get_primary_id(migrate_engine, table, cmd)
            if taxonomy_id is None:
                # Add a row to the taxonomy table.
                cmd = "INSERT INTO taxonomy VALUES (%s, %s, %s, '%s', '%s')"
                cmd = cmd % (nextval(migrate_engine, table),
                             localtimestamp(migrate_engine),
                             localtimestamp(migrate_engine),
                             species_name,
                             genus_name)
                migrate_engine.execute(cmd)
                taxonomy_table_inserts += 1
                taxonomy_id = get_latest_id(migrate_engine, table)

            # Process the experiment items.  Dependent tables: sample.
            table = "experiment"
            # See if we need to add a row to the experiment table.
            cmd = "SELECT id FROM experiment WHERE seq_facility = '%s' AND array_version = '%s'"
            cmd += " AND result_folder_name = '%s' and plate_barcode = '%s'"
            cmd = cmd % (seq_facility, array_version, result_folder_name, plate_barcode)
            experiment_id = get_primary_id(migrate_engine, table, cmd)
            if experiment_id is None:
                # Add a row to the experiment table.
                cmd = "INSERT INTO experiment VALUES (%s, %s, %s, '%s', '%s', '%s', '%s')"
                cmd = cmd % (nextval(migrate_engine, table),
                             localtimestamp(migrate_engine),
                             localtimestamp(migrate_engine),
                             seq_facility,
                             array_version,
                             result_folder_name,
                             plate_barcode)
                migrate_engine.execute(cmd)
                experiment_table_inserts += 1
                experiment_id = get_latest_id(migrate_engine, table)

            # Process the genotype items.  Dependent tables: sample.
            table = "genotype"
            # See if we need to add a row to the table.
            # Values for the following are not in the seed data.
            cmd = "SELECT id FROM genotype WHERE coral_mlg_clonal_id = '%s' AND coral_mlg_rep_sample_id = '%s'"
            cmd += " AND genetic_coral_species_call = '%s' AND bcoral_genet_id = '%s'"
            cmd = cmd % (coral_mlg_clonal_id,
                         coral_mlg_rep_sample_id,
                         genetic_coral_species_call,
                         bcoral_genet_id)
            genotype_id = get_primary_id(migrate_engine, table, cmd)
            if genotype_id is None:
                # Add a row to the table.
                cmd = "INSERT INTO genotype VALUES (%s, %s, %s, '%s', '%s', '%s', '%s')"
                cmd = cmd % (nextval(migrate_engine, table),
                             localtimestamp(migrate_engine),
                             localtimestamp(migrate_engine),
                             coral_mlg_clonal_id,
                             coral_mlg_rep_sample_id,
                             genetic_coral_species_call,
                             bcoral_genet_id)
                migrate_engine.execute(cmd)
                genotype_table_inserts += 1
                genotype_id = get_latest_id(migrate_engine, table)

            # Process the phenotype items.  Dependent tables: sample.
            table = "phenotype"
            # See if we need to add a row to the table.
            cmd = "SELECT id FROM phenotype WHERE disease_resist = '%s' AND bleach_resist = '%s' AND mortality = '%s'"
            cmd += " AND tle = '%s' AND spawning = '%s' AND sperm_motility = %s and healing_time = %s"
            cmd = cmd % (disease_resist,
                         bleach_resist,
                         mortality,
                         tle,
                         spawning,
                         sperm_motility,
                         healing_time)
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
                             sperm_motility,
                             healing_time)
                migrate_engine.execute(cmd)
                phenotype_table_inserts += 1
                phenotype_id = get_latest_id(migrate_engine, table)

            # Process the person items.  Dependent tables: sample.
            table = "person"
            # See if we need to add a row to the table.
            cmd = "SELECT id FROM person WHERE last_name = '%s' AND first_name = '%s' AND email = '%s'" % (collector_last_name, collector_first_name, email)
            person_id = get_primary_id(migrate_engine, table, cmd)
            if person_id is None:
                # Add a row to the table.
                cmd = "INSERT INTO person VALUES (%s, %s, %s, '%s', '%s', '%s', '%s')"
                cmd = cmd % (nextval(migrate_engine, table),
                             localtimestamp(migrate_engine),
                             localtimestamp(migrate_engine),
                             collector_last_name,
                             collector_first_name,
                             organization,
                             email)
                migrate_engine.execute(cmd)
                person_table_inserts += 1
                person_id = get_latest_id(migrate_engine, table)

            # Process the reef items.  Dependent tables: colony.
            table = "reef"
            # See if we need to add a row to the table.
            cmd = "SELECT id FROM reef WHERE name = '%s'" % reef
            reef_id = get_primary_id(migrate_engine, table, cmd)
            if reef_id is None:
                if geographic_origin == "reef":
                    # The geographic_origin value is used for deciding into which table
                    # to insert the latitude and longitude values.  If the geographic_origin
                    # is "reef", the values will be inserted into the reef table.
                    lat_val = latitude
                    long_val = longitude
                else:
                    lat_val = DEFAULT_MISSING_NUMERIC_VALUE
                    long_val = DEFAULT_MISSING_NUMERIC_VALUE
                # Add a row to the table.
                cmd = "INSERT INTO reef VALUES (%s, %s, %s, '%s', '%s', %s, %s, '%s')"
                cmd = cmd % (nextval(migrate_engine, table),
                             localtimestamp(migrate_engine),
                             localtimestamp(migrate_engine),
                             reef,
                             region,
                             lat_val,
                             long_val,
                             geographic_origin)
                migrate_engine.execute(cmd)
                reef_table_inserts += 1
                reef_id = get_latest_id(migrate_engine, table)

            # Process the colony items.  Dependent tables: sample.
            table = "colony"
            # See if we need to add a row to the table.
            cmd = "SELECT id FROM colony WHERE latitude = %s AND longitude = %s and reef_id = %s"
            cmd = cmd % (latitude, longitude, reef_id)
            colony_id = get_primary_id(migrate_engine, table, cmd)
            if colony_id is None:
                if geographic_origin == "colony":
                    # The geographic_origin value is used for deciding into which table
                    # to insert the latitude and longitude values.  If the geographic_origin
                    # is "colony", the values will be inserted into the colony table.
                    lat_val = latitude
                    long_val = longitude
                else:
                    lat_val = DEFAULT_MISSING_NUMERIC_VALUE
                    long_val = DEFAULT_MISSING_NUMERIC_VALUE
                # Add a row to the table.
                cmd = "INSERT INTO colony VALUES (%s, %s, %s, %s, %s, %s, %s)"
                cmd = cmd % (nextval(migrate_engine, table),
                             localtimestamp(migrate_engine),
                             localtimestamp(migrate_engine),
                             lat_val,
                             long_val,
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
                # We set allele_id top NULL here since we process the
                # allele data after the general seed data.  When the
                # allele data is processed, this column value will ne
                # updated with the foreign key value.
                allele_id = sql.null()
                percent_missing_data_sym = DEFAULT_MISSING_NUMERIC_VALUE
                percent_reference_sym = DEFAULT_MISSING_NUMERIC_VALUE
                percent_alternative_coral = DEFAULT_MISSING_NUMERIC_VALUE
                percent_heterozygous_sym = DEFAULT_MISSING_NUMERIC_VALUE
                # id, create_time, update_time, affy_id, sample_id,
                # allele_id, genotype_id, phenotype_id, experiment_id, colony_id,
                # colony_location, taxonomy_id, person_id
                cmd = "INSERT INTO sample VALUES (%s, %s, %s, '%s', '%s', %s, %s, %s, %s, %s, '%s', %s, %s, "
                if collection_date == "LOCALTIMESTAMP":
                    # collection_date
                    cmd += "%s, "
                else:
                    # collection_date
                    cmd += "'%s', "
                # user_specimen_id, registry_id, depth, dna_extraction_method, dna_concentration,
                # public, public_after_date, percent_missing_data_coral, percent_missing_data_sym, percent_reference_coral,
                # percent_reference_sym, percent_alternative_coral, percent_alternative_sym, percent_heterozygous_coral, percent_heterozygous_sym,
                # field_call
                cmd += "'%s', '%s', %s, '%s', %s, %s, '%s', %s, '%s', %s, '%s', '%s', '%s', %s, '%s', '%s')"
                cmd = cmd % (nextval(migrate_engine, table),
                             localtimestamp(migrate_engine),
                             localtimestamp(migrate_engine),
                             affy_id,
                             sample_id,
                             allele_id,
                             genotype_id,
                             phenotype_id,
                             experiment_id,
                             colony_id,
                             colony_location,
                             taxonomy_id,
                             person_id,
                             collection_date,
                             user_specimen_id,
                             registry_id,
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
                             percent_heterozygous_sym,
                             field_call)
                migrate_engine.execute(cmd)
                sample_table_inserts += 1
                sample_id = get_latest_id(migrate_engine, table)

    print("Inserted %d rows into the colony table." % colony_table_inserts)
    print("Inserted %d rows into the experiment table." % experiment_table_inserts)
    print("Inserted %d rows into the genotype table." % genotype_table_inserts)
    print("Inserted %d rows into the person table." % person_table_inserts)
    print("Inserted %d rows into the phenotype table." % phenotype_table_inserts)
    print("Inserted %d rows into the reef table." % reef_table_inserts)
    print("Inserted %d rows into the sample table." % sample_table_inserts)
    print("Inserted %d rows into the taxonomy table." % taxonomy_table_inserts)


def downgrade(migrate_engine):
    pass


def upgrade(migrate_engine):
    print(__doc__)
    metadata.bind = migrate_engine
    metadata.create_all()
    load_general_seed_data(migrate_engine)
    load_alleles_seed_data(migrate_engine)
    load_probe_annotation_table(migrate_engine)

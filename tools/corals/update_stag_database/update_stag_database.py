#!/usr/bin/env python
from __future__ import print_function

import argparse
import datetime
import dateutil.parser
import os
import psycopg2
import string
import subprocess
import sys

from six.moves import configparser

from sqlalchemy import create_engine
from sqlalchemy import MetaData
from sqlalchemy.engine.url import make_url

now = datetime.datetime.utcnow
metadata = MetaData()

DEFAULT_MISSING_NUMERIC_VALUE = -9.000000


def check_execution_errors(rc, fstderr, fstdout):
    if rc != 0:
        fh = open(fstdout, 'rb')
        out_msg = fh.read()
        fh.close()
        fh = open(fstderr, 'rb')
        err_msg = fh.read()
        fh.close()
        msg = '%s\n%s\n' % (str(out_msg), str(err_msg))
        sys.exit(msg)


def get_config_settings(config_file, section='defaults'):
    # Return a dictionary consisting of the key / value pairs
    # of the defaults section of config_file.
    d = {}
    config_parser = configparser.ConfigParser()
    config_parser.read(config_file)
    for key, value in config_parser.items(section):
        if section == 'defaults':
            d[string.upper(key)] = value
        else:
            d[key] = value
    return d


def get_response_buffers():
    fstderr = os.path.join(os.getcwd(), 'stderr.txt')
    fherr = open(fstderr, 'wb')
    fstdout = os.path.join(os.getcwd(), 'stdout.txt')
    fhout = open(fstdout, 'wb')
    return fstderr, fherr, fstdout, fhout


def get_sql_param_val_str(column_val, default):
    if set_to_null(column_val):
        val = default
    else:
        val = column_val
    return "= '%s'" % val


def get_value_from_config(config_defaults, value):
    return config_defaults.get(value, None)


def get_year_from_now():
    # Get current date plus one year for possible insertion
    # into the public_after_date column of the sample table.
    # The default behavior is for the value of the public
    # column to be True and the public_after_date to be NULL,
    # making the sample "public".  However, the user can
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
        return today.replace(year=year)
    except Exception:
        # Handle leap years.
        return today + (datetime.date(today.year + 1, 1, 1) - datetime.date(today.year, 1, 1))


def handle_column_value(val, get_sql_param=True, default=''):
    # Regarding the default value, a NULL value indicates an unknown value
    # and typically should not be confused with an empty string. Our application
    # does not need the concept of unknown value, so most columns are
    # non-nullable and our default is an empty string.
    param = handle_null(val)
    if get_sql_param:
        param_val_str = get_sql_param_val_str(val, default)
    if param is None:
        if get_sql_param:
            return default, param_val_str
        return default
    if get_sql_param:
        return param, param_val_str
    return param


def handle_null(val):
    if set_to_null(val):
        return None
    return val


def run_command(cmd):
    fstderr, fherr, fstdout, fhout = get_response_buffers()
    proc = subprocess.Popen(args=cmd, stderr=fherr, stdout=fhout, shell=True)
    rc = proc.wait()
    # Check results.
    fherr.close()
    fhout.close()
    check_execution_errors(rc, fstderr, fstdout)


def set_to_null(val):
    if val in ["", "NA", "NULL"]:
        return True
    return False


def split_line(line, sep="\t"):
    # Remove R quote chars.
    items = line.split(sep)
    unquoted_items = []
    for item in items:
        unquoted_items.append(item.strip('"'))
    return unquoted_items


def string_as_bool(string):
    if str(string).lower() in ('true', 'yes', 'on', '1'):
        return True
    else:
        return False


class StagDatabaseUpdater(object):
    def __init__(self):
        self.args = None
        self.conn = None
        self.parse_args()
        self.year_from_now = get_year_from_now()
        self.db_name = None
        self.db_storage_dir = None
        self.get_config_settings()
        self.outfh = open(self.args.output, "w")
        self.connect_db()
        self.engine = create_engine(self.args.database_connection_string)
        self.metadata = MetaData(self.engine)
        self.affy_ids = []
        self.allele_ids = []
        self.colony_ids = []
        self.experiment_ids = []
        self.genotype_ids = []
        self.person_ids = []
        self.phenotype_ids = []
        self.reef_ids = []
        self.taxonomy_ids = []

    def connect_db(self):
        url = make_url(self.args.database_connection_string)
        self.log('Connecting to database with URL: %s' % url)
        args = url.translate_connect_args(username='user')
        args.update(url.query)
        assert url.get_dialect().name == 'postgresql', 'This script can only be used with PostgreSQL.'
        self.conn = psycopg2.connect(**args)
        self.log("Successfully connected to the stag database...")

    def convert_date_string_for_database(self, date_string):
        # The value of date_string is %y/%m/%d with
        # the year being 2 digits (yikes!).
        fixed_century = "20%s" % date_string
        fixed_date = fixed_century.replace("/", "-")
        # Convert the string to a format required for
        # inserting into the database.
        database_format = dateutil.parser.parse(fixed_date)
        return str(database_format)

    def flush(self):
        self.conn.commit()

    def export_database(self):
        # Export the database to the configured storage location.
        if not os.path.isdir(self.db_storage_dir):
            os.makedirs(self.db_storage_dir)
        db_storage_path = os.path.join(self.db_storage_dir, "exported_%s_db" % self.db_name)
        cmd = "pg_dump %s -f %s" % (self.db_name, db_storage_path)
        run_command(cmd)

    def get_config_settings(self):
        config_defaults = get_config_settings(self.args.config_file)
        self.db_name = get_value_from_config(config_defaults, 'DB_NAME')
        self.db_storage_dir = get_value_from_config(config_defaults, 'DB_STORAGE_DIR')

    def get_next_sample_id(self):
        cmd = "SELECT sample_id FROM sample ORDER by id DESC;"
        cur = self.conn.cursor()
        cur.execute(cmd)
        try:
            last_sample_id = cur.fetchone()[0]
            # The value of last_sample_id will be something like A10171.
            last_sample_id_num = int(last_sample_id.lstrip("A"))
            next_sample_id_num = last_sample_id_num + 1
            next_sample_id = "A%d" % next_sample_id_num
        except Exception:
            next_sample_id = "A10000"
        return next_sample_id

    def log(self, msg):
        self.outfh.write("%s\n" % msg)

    def update_allele_table(self, file_path):
        self.log("Updating the allele table...")
        # Columns in the experiment file are:
        # affy_id allele
        allele_table_inserts = 0
        # The allele.tabular file contains a subset of the number of samples
        # to be inserted.  This is because those samples that failed will not
        # be included in the file.  Failed samples will have an affy_id value
        # of NA in self.affy_ids, which was generated when the genotype.tabular
        # file was processed, so we'll use that list to build the correct list
        # of self.allele_ids for later use when inserting into the sample table.
        fh = open(file_path, "r")
        # Skip the header
        fh.readline()
        for id_index, affy_id in enumerate(self.affy_ids):
            if set_to_null(affy_id):
                # This is a failed sample, so no allele strings will be
                # inserted, and we'll set the allele_id to the default
                # empty string.
                self.allele_ids.append("")
                continue
            # See if we need to add a row to the table.  The affy_id value
            # should not exist in the sample table.
            cmd = "SELECT allele_id FROM sample WHERE affy_id = '%s';" % affy_id
            cur = self.conn.cursor()
            cur.execute(cmd)
            try:
                allele_id = cur.fetchone()[0]
            except Exception:
                # Insert a row into the allele table.
                line = fh.readline()
                line = line.rstrip()
                items = split_line(line)
                allele = items[1]
                cmd = "INSERT INTO allele VALUES (nextval('allele_id_seq'), %s, %s, %s) RETURNING id;"
                args = ['NOW()', 'NOW()', allele]
                cur = self.update(cmd, args)
                self.flush()
                allele_id = cur.fetchone()[0]
                allele_table_inserts += 1
            self.allele_ids.append(allele_id)
        self.log("Inserted %d rows into the allele table..." % allele_table_inserts)

    def update_colony_table(self, file_path):
        self.log("Updating the colony table...")
        # Columns in the colony file are:
        # latitude longitude depth geographic_origin
        # The geographic_origin value is used for deciding into which table
        # to insert the latitude and longitude values.  If the geographic_origin
        # is "colony", the values will be inserted into the colony table.
        colony_table_inserts = 0
        with open(file_path) as fh:
            for i, line in enumerate(fh):
                if i == 0:
                    # Skip header
                    continue
                # Keep track of foreign keys since we skip the header line.
                id_index = i - 1
                line = line.rstrip()
                items = split_line(line)
                geographic_origin = items[3]
                if set_to_null(geographic_origin):
                    geographic_origin = "reef"
                else:
                    geographic_origin = geographic_origin.lower()
                if geographic_origin == "colony":
                    latitude = "%6f" % float(items[0])
                    longitude = "%6f" % float(items[1])
                else:
                    latitude = DEFAULT_MISSING_NUMERIC_VALUE
                    longitude = DEFAULT_MISSING_NUMERIC_VALUE
                depth = handle_column_value(items[2], get_sql_param=False, default=-9)
                reef_id = self.reef_ids[id_index]
                # See if we need to add a row to the table.
                cmd = "SELECT id FROM colony WHERE latitude = %s " % latitude
                cmd += "AND longitude = %s AND depth = %s " % (longitude, depth)
                cmd += "AND reef_id = %s;" % reef_id
                cur = self.conn.cursor()
                cur.execute(cmd)
                try:
                    colony_id = cur.fetchone()[0]
                except Exception:
                    # Insert a row into the colony table.
                    cmd = "INSERT INTO colony VALUES (nextval('colony_id_seq'), %s, %s, %s, %s, %s, %s) RETURNING id;"
                    args = ['NOW()', 'NOW()', longitude, latitude, depth, reef_id]
                    cur = self.update(cmd, args)
                    self.flush()
                    colony_id = cur.fetchone()[0]
                    colony_table_inserts += 1
                self.colony_ids.append(colony_id)
        self.log("Inserted %d rows into the colony table..." % colony_table_inserts)

    def update_experiment_table(self, file_path):
        self.log("Updating the experiment table...")
        # Columns in the experiment file are:
        # seq_facility array_version result_folder_name plate_barcode
        experiment_table_inserts = 0
        with open(file_path) as fh:
            for i, line in enumerate(fh):
                if i == 0:
                    # Skip header
                    continue
                # Keep track of foreign keys since we skip the header line.
                line = line.rstrip()
                items = split_line(line)
                seq_facility, seq_facility_param_val_str = handle_column_value(items[0])
                array_version, array_version_param_val_str = handle_column_value(items[1])
                result_folder_name, result_folder_name_param_val_str = handle_column_value(items[2])
                plate_barcode, plate_barcode_param_val_str = handle_column_value(items[3])
                # See if we need to add a row to the table.
                cmd = "SELECT id FROM experiment WHERE seq_facility %s " % seq_facility_param_val_str
                cmd += "AND array_version %s " % array_version_param_val_str
                cmd += "AND result_folder_name %s " % result_folder_name_param_val_str
                cmd += "AND plate_barcode %s;" % plate_barcode_param_val_str
                cur = self.conn.cursor()
                cur.execute(cmd)
                try:
                    experiment_id = cur.fetchone()[0]
                except Exception:
                    # Insert a row into the experiment table.
                    cmd = "INSERT INTO experiment VALUES (nextval('experiment_id_seq'), %s, %s, %s, %s, %s, %s) RETURNING id;"
                    args = ['NOW()', 'NOW()', seq_facility, array_version, result_folder_name, plate_barcode]
                    cur = self.update(cmd, args)
                    self.flush()
                    experiment_id = cur.fetchone()[0]
                    experiment_table_inserts += 1
                self.experiment_ids.append(experiment_id)
        self.log("Inserted %d rows into the experiment table..." % experiment_table_inserts)

    def update_genotype_table(self, file_path):
        self.log("Updating the genotype table...")
        # Columns in the genotype file are:
        # affy_id coral_mlg_clonal_id user_specimen_id db_match genetic_coral_species_call
        # coral_mlg_rep_sample_id
        genotype_table_inserts = 0
        with open(file_path) as fh:
            for i, line in enumerate(fh):
                if i == 0:
                    # Skip header
                    continue
                line = line.rstrip()
                items = split_line(line)
                # Keep an in-memory list of affy_ids for use
                # when updating the allele table.
                self.affy_ids.append(items[0])
                coral_mlg_clonal_id = items[1]
                # The value of db_match will be "no_match" if
                # a new row should be inserted into the table.
                db_match = items[3].lower()
                genetic_coral_species_call = handle_column_value(items[4], get_sql_param=False)
                coral_mlg_rep_sample_id = handle_column_value(items[5], get_sql_param=False)
                if db_match == "failed":
                    # Handle the special case of a failed sample.
                    cmd = "SELECT id FROM genotype WHERE coral_mlg_clonal_id = 'failed'"
                else:
                    cmd = "SELECT id FROM genotype WHERE coral_mlg_clonal_id = '%s'" % coral_mlg_clonal_id
                # See if we need to add a row to the table.
                cur = self.conn.cursor()
                cur.execute(cmd)
                try:
                    genotype_id = cur.fetchone()[0]
                    if db_match == "failed":
                        val = db_match
                    else:
                        val = "match"
                    self.log("Found genotype row with id %d, value of db_match: %s, should be %s." % (genotype_id, db_match, val))
                except Exception:
                    # Insert a row into the genotype table.
                    cmd = "INSERT INTO genotype VALUES (nextval('genotype_id_seq'), NOW(), NOW(), "
                    cmd += "'%s', '%s', '%s') RETURNING id;"
                    cmd = cmd % (coral_mlg_clonal_id, coral_mlg_rep_sample_id, genetic_coral_species_call)
                    args = [coral_mlg_clonal_id, coral_mlg_rep_sample_id, genetic_coral_species_call]
                    cur = self.update(cmd, args)
                    self.flush()
                    genotype_id = cur.fetchone()[0]
                    if db_match == "failed":
                        val = db_match
                    else:
                        val = "no_match"
                    self.log("Inserted genotype row with id %d, value of db_match: %s, should be %s." % (genotype_id, db_match, val))
                    genotype_table_inserts += 1
                self.genotype_ids.append(genotype_id)
        self.log("Inserted %d rows into the genotype table..." % genotype_table_inserts)

    def update_person_table(self, file_path):
        self.log("Updating the person table...")
        # Columns in the person file are:
        # last_name first_name organization email
        person_table_inserts = 0
        with open(file_path) as fh:
            for i, line in enumerate(fh):
                if i == 0:
                    # Skip header
                    continue
                line = line.rstrip()
                items = split_line(line)
                last_name = items[0]
                first_name = items[1]
                organization = items[2]
                email = items[3]
                # See if we need to add a row to the table.
                cmd = "SELECT id FROM person WHERE email = '%s';" % email
                cur = self.conn.cursor()
                cur.execute(cmd)
                try:
                    person_id = cur.fetchone()[0]
                except Exception:
                    # Insert a row into the person table.
                    cmd = "INSERT INTO person VALUES (nextval('person_id_seq'), NOW(), NOW(), "
                    cmd += "%s, %s, %s, %s) RETURNING id;"
                    args = [last_name, first_name, organization, email]
                    cur = self.update(cmd, args)
                    self.flush()
                    person_id = cur.fetchone()[0]
                    person_table_inserts += 1
                self.person_ids.append(person_id)
        self.log("Inserted %d rows into the person table..." % person_table_inserts)

    def update_phenotype_table(self, file_path):
        self.log("Updating the phenotype table...")
        # Columns in the phenotype file are:
        # disease_resist bleach_resist mortality tle spawning sperm_motility healing_time
        phenotype_table_inserts = 0
        with open(file_path) as fh:
            for i, line in enumerate(fh):
                if i == 0:
                    # Skip header
                    continue
                line = line.rstrip()
                items = split_line(line)
                disease_resist, disease_resist_param_val_str = handle_column_value(items[0], default=-9)
                bleach_resist, bleach_resist_param_val_str = handle_column_value(items[1], default=-9)
                mortality, mortality_param_val_str = handle_column_value(items[2], default=-9)
                tle, tle_param_val_str = handle_column_value(items[3], default=-9)
                spawning, spawning_param_val_str = handle_column_value(items[4])
                sperm_motility, sperm_motility_param_val_str = handle_column_value(items[5], default=-9.0)
                healing_time, healing_time_param_val_str = handle_column_value(items[6], default=-9.0)
                # See if we need to add a row to the phenotype table.
                cmd = " SELECT id FROM phenotype WHERE disease_resist %s "
                cmd += "AND bleach_resist %s AND mortality %s AND tle %s "
                cmd += "AND spawning %s AND sperm_motility %s AND healing_time %s;"
                cmd = cmd % (disease_resist_param_val_str, bleach_resist_param_val_str,
                             mortality_param_val_str, tle_param_val_str, spawning_param_val_str,
                             sperm_motility_param_val_str, healing_time_param_val_str)
                cur = self.conn.cursor()
                cur.execute(cmd)
                try:
                    phenotype_id = cur.fetchone()[0]
                except Exception:
                    # Insert a row into the phenotype table.
                    cmd = "INSERT INTO phenotype VALUES (nextval('phenotype_id_seq'), NOW(), NOW(), "
                    cmd += "%s, %s, %s, %s, %s, %s, %s) RETURNING id;"
                    args = [disease_resist, bleach_resist, mortality, tle, spawning, sperm_motility, healing_time]
                    cur = self.update(cmd, args)
                    self.flush()
                    phenotype_id = cur.fetchone()[0]
                    phenotype_table_inserts += 1
                self.phenotype_ids.append(phenotype_id)
        self.log("Inserted %d rows into the phenotype table..." % phenotype_table_inserts)

    def update_reef_table(self, file_path):
        self.log("Updating the reef table...")
        # Columns in the reef file are:
        # name region latitude longitude geographic_origin
        # The geographic_origin value is used for deciding into which table
        # to insert the latitude and longitude values.  If the geographic_origin
        # is "reef", the values will be inserted into the reef table.
        reef_table_inserts = 0
        with open(file_path) as fh:
            for i, line in enumerate(fh):
                if i == 0:
                    # Skip header
                    continue
                line = line.rstrip()
                items = split_line(line)
                name = items[0]
                region = items[1]
                geographic_origin = items[4]
                if set_to_null(geographic_origin):
                    geographic_origin = "reef"
                else:
                    geographic_origin = geographic_origin.lower()
                if geographic_origin == "reef":
                    latitude = "%6f" % float(items[2])
                    longitude = "%6f" % float(items[3])
                else:
                    latitude = DEFAULT_MISSING_NUMERIC_VALUE
                    longitude = DEFAULT_MISSING_NUMERIC_VALUE
                # See if we need to add a row to the reef table.
                cmd = "SELECT id FROM reef WHERE name = $$%s$$ AND region = '%s' " % (name, region)
                cmd += "AND latitude = %s AND longitude = %s " % (latitude, longitude)
                cmd += "AND geographic_origin = '%s';" % geographic_origin
                cur = self.conn.cursor()
                cur.execute(cmd)
                try:
                    reef_id = cur.fetchone()[0]
                except Exception:
                    # Insert a row into the reef table.
                    cmd = "INSERT INTO reef VALUES (nextval('reef_id_seq'), %s, %s, %s, %s, %s, %s, %s) RETURNING id;"
                    args = ['NOW()', 'NOW()', name, region, latitude, longitude, geographic_origin]
                    cur = self.update(cmd, args)
                    self.flush()
                    reef_id = cur.fetchone()[0]
                    reef_table_inserts += 1
                self.reef_ids.append(reef_id)
        self.log("Inserted %d rows into the reef table..." % reef_table_inserts)

    def update_sample_table(self, file_path):
        self.log("Updating the sample table...")
        # Columns in the sample file are:
        # affy_id colony_location collection_date user_specimen_id registry_id
        # depth dna_extraction_method dna_concentration public public_after_date
        # percent_missing_data_coral percent_missing_data_sym percent_reference_coral percent_reference_sym percent_alternative_coral
        # percent_alternative_sym percent_heterozygous_coral percent_heterozygous_sym field_call, bcoral_genet_id
        sample_table_inserts = 0
        with open(file_path) as fh:
            for i, line in enumerate(fh):
                if i == 0:
                    # Skip header
                    continue
                line = line.rstrip()
                # Keep track of foreign keys since we skip the header line.
                id_index = i - 1
                items = split_line(line)
                sample_id = self.get_next_sample_id()
                allele_id = self.allele_ids[id_index]
                genotype_id = self.genotype_ids[id_index]
                phenotype_id = self.phenotype_ids[id_index]
                experiment_id = self.experiment_ids[id_index]
                colony_id = self.colony_ids[id_index]
                colony_location = handle_column_value(items[1], get_sql_param=False)
                taxonomy_id = self.taxonomy_ids[id_index]
                collector_id = self.person_ids[id_index]
                collection_date = items[2]
                user_specimen_id = items[3]
                affy_id = handle_column_value(items[0], get_sql_param=False, default="%s_%s" % (sample_id, user_specimen_id))
                registry_id = handle_column_value(items[4], get_sql_param=False, default=-9)
                depth = handle_column_value(items[5], get_sql_param=False, default=-9)
                dna_extraction_method = handle_column_value(items[6], get_sql_param=False)
                dna_concentration = handle_column_value(items[7], get_sql_param=False)
                public = items[8]
                if string_as_bool(public):
                    public_after_date = ''
                else:
                    if set_to_null(items[9]):
                        public_after_date = self.year_from_now
                    else:
                        public_after_date = items[9]
                percent_missing_data_coral = handle_column_value(items[10], get_sql_param=False)
                percent_missing_data_sym = handle_column_value(items[11], get_sql_param=False)
                percent_reference_coral = handle_column_value(items[12], get_sql_param=False)
                percent_reference_sym = handle_column_value(items[13], get_sql_param=False)
                percent_alternative_coral = handle_column_value(items[14], get_sql_param=False)
                percent_alternative_sym = handle_column_value(items[15], get_sql_param=False)
                percent_heterozygous_coral = handle_column_value(items[16], get_sql_param=False)
                percent_heterozygous_sym = handle_column_value(items[17], get_sql_param=False)
                field_call = handle_column_value(items[18], get_sql_param=False)
                bcoral_genet_id = handle_column_value(items[19], get_sql_param=False)
                # Insert a row into the sample table.
                cmd = "INSERT INTO sample VALUES (nextval('sample_id_seq'), %s, %s, %s, %s, %s, %s, %s, "
                cmd += "%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, "
                cmd += "%s, %s, %s, %s) RETURNING id;"
                args = ['NOW()', 'NOW()', affy_id, sample_id, allele_id, genotype_id, phenotype_id,
                        experiment_id, colony_id, colony_location, taxonomy_id, collector_id,
                        collection_date, user_specimen_id, registry_id, depth,
                        dna_extraction_method, dna_concentration, public, public_after_date,
                        percent_missing_data_coral, percent_missing_data_sym, percent_reference_coral,
                        percent_reference_sym, percent_alternative_coral, percent_alternative_sym,
                        percent_heterozygous_coral, percent_heterozygous_sym, field_call, bcoral_genet_id]
                cur = self.update(cmd, args)
                self.flush()
                sample_id = cur.fetchone()[0]
                sample_table_inserts += 1
        self.log("Inserted %d rows into the sample table..." % sample_table_inserts)

    def update_taxonomy_table(self, file_path):
        self.log("Updating the taxonomy table...")
        # Columns in the taxonomy file are:
        # genetic_coral_species_call affy_id genus_name species_name"
        taxonomy_table_inserts = 0
        with open(file_path) as fh:
            for i, line in enumerate(fh):
                if i == 0:
                    # Skip header
                    continue
                line = line.rstrip()
                items = split_line(line)
                genus_name = handle_column_value(items[2], get_sql_param=False, default='unknown')
                species_name = handle_column_value(items[3], get_sql_param=False, default='unknown')
                # See if we need to add a row to the taxonomy table.
                cmd = "SELECT id FROM taxonomy WHERE species_name = '%s' AND genus_name = '%s';" % (species_name, genus_name)
                cur = self.conn.cursor()
                cur.execute(cmd)
                try:
                    taxonomy_id = cur.fetchone()[0]
                except Exception:
                    # Insert a row into the taxonomy table.
                    cmd = "INSERT INTO taxonomy VALUES (nextval('taxonomy_id_seq'), %s, %s, %s, %s) RETURNING id;"
                    args = ['NOW()', 'NOW()', species_name, genus_name]
                    cur = self.update(cmd, args)
                    self.flush()
                    taxonomy_id = cur.fetchone()[0]
                    taxonomy_table_inserts += 1
                self.taxonomy_ids.append(taxonomy_id)
        self.log("Inserted %d rows into the taxonomy table..." % taxonomy_table_inserts)

    def parse_args(self):
        parser = argparse.ArgumentParser()
        parser.add_argument('--config_file', dest='config_file', help='usd_config.ini'),
        parser.add_argument('--database_connection_string', dest='database_connection_string', help='Postgres database connection string'),
        parser.add_argument('--input_dir', dest='input_dir', help='Input datasets for database insertion')
        parser.add_argument('--output', dest='output', help='Output dataset'),
        self.args = parser.parse_args()

    def run(self):
        self.export_database()
        input_dir = self.args.input_dir
        for file_name in os.listdir(input_dir):
            # Tables must be loaded in such a way that foreign keys
            # are properly handled.  The sample table must be loaded
            # last.
            if file_name.startswith("allele"):
                allele_file = os.path.join(input_dir, file_name)
            if file_name.startswith("colony"):
                colony_file = os.path.join(input_dir, file_name)
            if file_name.startswith("experiment"):
                experiment_file = os.path.join(input_dir, file_name)
            if file_name.startswith("genotype"):
                genotype_file = os.path.join(input_dir, file_name)
            elif file_name.startswith("person"):
                person_file = os.path.join(input_dir, file_name)
            elif file_name.startswith("phenotype"):
                phenotype_file = os.path.join(input_dir, file_name)
            elif file_name.startswith("reef"):
                reef_file = os.path.join(input_dir, file_name)
            elif file_name.startswith("sample"):
                sample_file = os.path.join(input_dir, file_name)
            elif file_name.startswith("taxonomy"):
                taxonomy_file = os.path.join(input_dir, file_name)
        # Now tables can be loaded in the appropriate order.
        self.update_experiment_table(experiment_file)
        self.update_genotype_table(genotype_file)
        self.update_allele_table(allele_file)
        self.update_person_table(person_file)
        self.update_phenotype_table(phenotype_file)
        self.update_reef_table(reef_file)
        self.update_colony_table(colony_file)
        self.update_taxonomy_table(taxonomy_file)
        self.update_sample_table(sample_file)

    def shutdown(self):
        self.log("Shutting down...")
        self.conn.close()

    def stop_err(self, msg):
        sys.stderr.write(msg)
        self.outfh.flush()
        self.outfh.close()
        sys.exit(1)

    def update(self, cmd, args):
        for i, arg in enumerate(args):
            args[i] = handle_null(arg)
        try:
            cur = self.conn.cursor()
            cur.execute(cmd, tuple(args))
        except Exception as e:
            msg = "Caught exception executing SQL:\n%s\nException:\n%s\n" % (cmd.format(args), e)
            self.stop_err(msg)
        return cur


if __name__ == '__main__':
    sdu = StagDatabaseUpdater()
    sdu.run()
    sdu.shutdown()

#!/usr/bin/env python


import argparse
import datetime
import os
import psycopg2
import subprocess
import sys

from six.moves import configparser

from sqlalchemy import create_engine
from sqlalchemy import MetaData
from sqlalchemy.engine.url import make_url

now = datetime.datetime.utcnow
metadata = MetaData()


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
            d[key.upper()] = value
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


class StagDatabaseUpdater(object):
    def __init__(self):
        self.args = None
        self.conn = None
        self.parse_args()
        self.db_name = None
        self.db_storage_dir = None
        self.get_config_settings()
        self.outfh = open(self.args.output, "w")
        self.connect_db()
        self.engine = create_engine(self.args.database_connection_string)
        self.metadata = MetaData(self.engine)
        self.affy_ids = []

    def connect_db(self):
        url = make_url(self.args.database_connection_string)
        self.log('Attempting to connect to the database...\n')
        args = url.translate_connect_args(username='user')
        args.update(url.query)
        assert url.get_dialect().name == 'postgresql', 'This script can only be used with PostgreSQL.'
        self.conn = psycopg2.connect(**args)
        self.log("Successfully connected to the database...\n")

    def export_database(self):
        # Export the database to the configured storage location.
        if not os.path.isdir(self.db_storage_dir):
            os.makedirs(self.db_storage_dir)
        db_storage_path = os.path.join(self.db_storage_dir, "exported_%s_db" % self.db_name)
        cmd = "pg_dump %s -f %s" % (self.db_name, db_storage_path)
        run_command(cmd)

    def flush(self):
        self.conn.commit()

    def get_config_settings(self):
        config_defaults = get_config_settings(self.args.config_file)
        self.db_name = get_value_from_config(config_defaults, 'DB_NAME')
        base_storage_dir = get_value_from_config(config_defaults, 'DB_STORAGE_DIR')
        # Use the date to name the storage directory to
        # enable storing a file per day (multiple runs
        # per day will overwrite the existing file.
        date_str = datetime.datetime.now().strftime("%Y_%m_%d")
        self.db_storage_dir = os.path.join(base_storage_dir, date_str)

    def handle_none(self, val):
        if val.lower() in ['none']:
            return ''
        return val

    def log(self, msg):
        self.outfh.write("%s\n" % msg)

    def parse_args(self):
        parser = argparse.ArgumentParser()
        parser.add_argument('--config_file', dest='config_file', help='usd_config.ini'),
        parser.add_argument('--database_connection_string', dest='database_connection_string', help='Postgres database connection string'),
        parser.add_argument('--input', dest='input', help='Input tabular file containing entries for updating the stag databse')
        parser.add_argument('--output', dest='output', help='Output dataset'),
        self.args = parser.parse_args()

    def parse_line(self, i, line):
        # Columns in input
        # [0]affy_id [1]sample_id [2]user_specimen_id [3]field_call [4]depth
        # [5]percent_missing_data_coral [6]percent_heterozygous_coral [7]percent_acerv_coral [8]percent_apalm_coral [9]bcoral_genet_id
        # [10]registry_id [11]dna_extraction_method [12]dna_concentration [13]colony_location [14]colony_latitude
        # [15]colony_longitude [16]colony_depth [17]reef_name [18]region [19]reef_latitude
        # [20]reef_longitude [21]geographic_origin [22]coral_mlg_clonal_id [23]coral_mlg_rep_sample_id
        # [24]genetic_coral_species_call [25]spawning [26]sperm_motility [27]tle [28]disease_resist [29]bleach_resist
        # [30]mortality [31]healing_time [32]seq_facility [33]array_version [34]plate_barcode
        # [35]collector_last_name [36]collector_first_name [37]organization [38]email [39]collection_date
        items = line.split("\t")
        if len(items) != 40:
            # Skip bad lines.
            self.outfh.write("\nSkipping line %d, shown below, with %d items, 40 items are required.\n%s\n\n" % (i + 1, len(items), line))
            return (), (), (), (), (), (), ()
        affy_id = self.handle_none(items[0])
        # The sample_id column value is i
        # auto-generated so should not be
        # changed in this tool..
        # sample_id = self.handle_none(items[1])
        user_specimen_id = self.handle_none(items[2])
        field_call = self.handle_none(items[3])
        depth = self.handle_none(items[4])
        percent_missing_data_coral = self.handle_none(items[5])
        percent_heterozygous_coral = self.handle_none(items[6])
        percent_acerv_coral = self.handle_none(items[7])
        percent_apalm_coral = self.handle_none(items[8])
        bcoral_genet_id = self.handle_none(items[9])
        registry_id = self.handle_none(items[10])
        dna_extraction_method = self.handle_none(items[11])
        dna_concentration = self.handle_none(items[12])
        colony_location = self.handle_none(items[13])
        colony_latitude = self.handle_none(items[14])
        colony_longitude = self.handle_none(items[15])
        colony_depth = self.handle_none(items[16])
        reef_name = self.handle_none(items[17])
        region = self.handle_none(items[18])
        reef_latitude = self.handle_none(items[19])
        reef_longitude = self.handle_none(items[20])
        geographic_origin = self.handle_none(items[21])
        if len(geographic_origin) > 0:
            geographic_origin = geographic_origin.lower()
        coral_mlg_clonal_id = self.handle_none(items[22])
        coral_mlg_rep_sample_id = self.handle_none(items[23])
        genetic_coral_species_call = self.handle_none(items[24])
        spawning = self.handle_none(items[25])
        sperm_motility = self.handle_none(items[26])
        tle = self.handle_none(items[27])
        disease_resist = self.handle_none(items[28])
        bleach_resist = self.handle_none(items[29])
        mortality = self.handle_none(items[30])
        healing_time = self.handle_none(items[31])
        seq_facility = self.handle_none(items[32])
        array_version = self.handle_none(items[33])
        plate_barcode = self.handle_none(items[34])
        collector_last_name = self.handle_none(items[35])
        collector_first_name = self.handle_none(items[36])
        organization = self.handle_none(items[37])
        email = self.handle_none(items[38])
        if len(email) > 0:
            email = email.lower()
        collection_date = self.handle_none(items[39])
        if not self.valid_date(collection_date):
            self.log("\nSkipping line %d, shown below, due to invalid collection_date, date formats must be yyyy-mm-dd.\n%s\n\n" % (i + 1, line))
            return (), (), (), (), (), (), ()
        # Gather values per table columns.
        colony_tup = (colony_latitude, colony_longitude, colony_depth)
        experiment_tup = (seq_facility, array_version, plate_barcode)
        genotype_tup = (coral_mlg_clonal_id, coral_mlg_rep_sample_id, genetic_coral_species_call)
        person_tup = (collector_last_name, collector_first_name, organization, email)
        phenotype_tup = (spawning, sperm_motility, tle, disease_resist, bleach_resist, mortality, healing_time)
        reef_tup = (reef_name, region, reef_latitude, reef_longitude, geographic_origin)
        sample_tup = (affy_id, user_specimen_id, field_call, depth, percent_missing_data_coral,
                      percent_heterozygous_coral, percent_acerv_coral, percent_apalm_coral, bcoral_genet_id,
                      registry_id, dna_extraction_method, dna_concentration, colony_location, collection_date)
        return colony_tup, experiment_tup, genotype_tup, person_tup, phenotype_tup, reef_tup, sample_tup

    def run(self):
        self.export_database()
        with open(self.args.input, "r") as fh:
            for i, line in enumerate(fh):
                line = line.rstrip()
                if i == 0:
                    # Skip header.
                    continue
                colony_tup, experiment_tup, genotype_tup, person_tup, phenotype_tup, reef_tup, sample_tup = self.parse_line(i, line)
                if len(sample_tup) == 0:
                    # Skip invalid lines.
                    continue
                # The affy_id is the unique key in the sample
                # table that enables connections to other tables.
                affy_id = sample_tup[0]
                self.log("\nProcessing affy id: %s...\n" % affy_id)

                # Get the foreign key ids from the sample table.
                cmd = "SELECT colony_id, experiment_id, genotype_id, collector_id, phenotype_id "
                cmd += "FROM sample WHERE affy_id = '%s';" % affy_id
                cur = self.conn.cursor()
                cur.execute(cmd)
                try:
                    vals = cur.fetchone()
                    colony_id = vals[0]
                    experiment_id = vals[1]
                    genotype_id = vals[2]
                    collector_id = vals[3]
                    phenotype_id = vals[4]
                except Exception as e:
                    msg = "Caught exception executing SQL:\n%s\nException:\n%s\n" % (cmd, e)
                    sys.stderr.write(msg)
                    self.outfh.flush()
                    self.outfh.close()
                    self.conn.close()
                    sys.exit(1)
                self.update_colony_table(colony_id, affy_id, colony_tup)
                self.update_experiment_table(experiment_id, affy_id, experiment_tup)
                self.update_genotype_table(genotype_id, affy_id, genotype_tup)
                self.update_person_table(collector_id, affy_id, person_tup)
                self.update_phenotype_table(phenotype_id, affy_id, phenotype_tup)
                self.update_reef_table(colony_id, affy_id, reef_tup)
                self.update_sample_table(sample_tup)

    def shutdown(self):
        self.log("\nShutting down...\n")
        self.outfh.flush()
        self.outfh.close()
        self.conn.close()

    def update(self, cmd, args):
        for i, arg in enumerate(args):
            args[i] = handle_null(arg)
        try:
            cur = self.conn.cursor()
            cur.execute(cmd, tuple(args))
        except Exception as e:
            msg = "Caught exception executing SQL:\n%s\nException:\n%s\n" % (cmd.format(args), e)
            sys.stderr.write(msg)
            self.outfh.flush()
            self.outfh.close()
            self.conn.close()
            sys.exit(1)

    def update_colony_table(self, colony_id, affy_id, tup):
        latitude = "%6f" % float(tup[0])
        longitude = "%6f" % float(tup[1])
        depth = handle_column_value(tup[2], get_sql_param=False, default=-9.0)
        # Update the row in the colony table.
        cmd = "UPDATE colony SET latitude=%s, longitude=%s, depth=%s WHERE id = %s;"
        args = [latitude, longitude, depth, colony_id]
        self.update(cmd, args)
        self.flush()
        self.log("Updated the colony table row associated with affy_id '%s'\n" % affy_id)
        return colony_id

    def update_experiment_table(self, experiment_id, affy_id, tup):
        seq_facility = handle_column_value(tup[0], get_sql_param=False)
        array_version = handle_column_value(tup[1], get_sql_param=False)
        plate_barcode = handle_column_value(tup[2], get_sql_param=False)
        # Update the row in the experiment table.
        cmd = "UPDATE experiment SET seq_facility=%s, array_version=%s, plate_barcode=%s WHERE id=%s;"
        args = [seq_facility, array_version, plate_barcode, experiment_id]
        self.update(cmd, args)
        self.flush()
        self.log("Updated the experiment table row associated with affy_id '%s'\n" % affy_id)

    def update_genotype_table(self, genotype_id, affy_id, tup):
        coral_mlg_clonal_id = handle_column_value(tup[0], get_sql_param=False)
        coral_mlg_rep_sample_id = handle_column_value(tup[1], get_sql_param=False)
        genetic_coral_species_call = handle_column_value(tup[2], get_sql_param=False)
        # Update the row in the genotype table.
        cmd = "UPDATE genotype SET coral_mlg_clonal_id=%s, coral_mlg_rep_sample_id=%s, genetic_coral_species_call=%s WHERE id=%s;"
        args = [coral_mlg_clonal_id, coral_mlg_rep_sample_id, genetic_coral_species_call, genotype_id]
        self.update(cmd, args)
        self.flush()
        self.log("Updated the genotype table row associated with affy_id '%s'\n" % affy_id)

    def update_person_table(self, collector_id, affy_id, tup):
        last_name = handle_column_value(tup[0], get_sql_param=False)
        first_name = handle_column_value(tup[1], get_sql_param=False)
        organization = handle_column_value(tup[2], get_sql_param=False)
        email = handle_column_value(tup[3], get_sql_param=False)
        # Update the row in the person table.
        cmd = "UPDATE person SET last_name=%s, first_name=%s, organization=%s, email=%s WHERE id=%s;"
        args = [last_name, first_name, organization, email, collector_id]
        self.update(cmd, args)
        self.flush()
        self.log("Updated the person table row associated with affy_id '%s'\n" % affy_id)

    def update_phenotype_table(self, phenotype_id, affy_id, tup):
        spawning = handle_column_value(tup[0], get_sql_param=False)
        sperm_motility = handle_column_value(tup[1], get_sql_param=False)
        tle = handle_column_value(tup[2], get_sql_param=False)
        disease_resist = handle_column_value(tup[3], get_sql_param=False)
        bleach_resist = handle_column_value(tup[4], get_sql_param=False)
        mortality = handle_column_value(tup[5], get_sql_param=False)
        healing_time = handle_column_value(tup[6], get_sql_param=False)
        # Update the row in the phenotype table.
        cmd = "UPDATE phenotype SET spawning=%s, sperm_motility=%s, tle=%s, disease_resist=%s, "
        cmd += "bleach_resist=%s, mortality=%s, healing_time=%s WHERE id=%s;"
        args = [spawning, sperm_motility, tle, disease_resist, bleach_resist, mortality, healing_time, phenotype_id]
        self.update(cmd, args)
        self.flush()
        self.log("Updated the phenotype table row associated with affy_id '%s'\n" % affy_id)

    def update_reef_table(self, colony_id, affy_id, tup):
        # Get the reef_id key ids from the colony table.
        cmd = "SELECT reef_id FROM colony WHERE id = %s;" % colony_id
        cur = self.conn.cursor()
        cur.execute(cmd)
        try:
            reef_id = cur.fetchone()[0]
        except Exception as e:
            msg = "Caught exception executing SQL:\n%s\nException:\n%s\n" % (cmd, e)
            sys.stderr.write(msg)
            self.outfh.flush()
            self.outfh.close()
            self.conn.close()
            sys.exit(1)
        name = handle_column_value(tup[0], get_sql_param=False)
        region = handle_column_value(tup[1], get_sql_param=False)
        latitude = handle_column_value(tup[2], get_sql_param=False)
        longitude = handle_column_value(tup[3], get_sql_param=False)
        geographic_origin = handle_column_value(tup[4], get_sql_param=False)
        # Update the row in the reef table.
        cmd = "UPDATE reef SET name=%s, region=%s, latitude=%s, longitude=%s, geographic_origin=%s WHERE id=%s;"
        args = [name, region, latitude, longitude, geographic_origin, reef_id]
        self.update(cmd, args)
        self.flush()
        self.log("Updated the reef table row associated with affy_id '%s'\n" % affy_id)

    def update_sample_table(self, tup):
        affy_id = tup[0]
        user_specimen_id = handle_column_value(tup[1], get_sql_param=False)
        field_call = handle_column_value(tup[2], get_sql_param=False)
        depth = handle_column_value(tup[3], get_sql_param=False, default=-9.0)
        percent_missing_data_coral = handle_column_value(tup[4], get_sql_param=False)
        percent_heterozygous_coral = handle_column_value(tup[5], get_sql_param=False)
        percent_acerv_coral = handle_column_value(tup[6], get_sql_param=False)
        percent_apalm_coral = handle_column_value(tup[7], get_sql_param=False)
        bcoral_genet_id = handle_column_value(tup[8], get_sql_param=False)
        registry_id = handle_column_value(tup[9], get_sql_param=False, default=-9)
        dna_extraction_method = handle_column_value(tup[10], get_sql_param=False)
        dna_concentration = handle_column_value(tup[11], get_sql_param=False)
        colony_location = handle_column_value(tup[12], get_sql_param=False)
        collection_date = tup[13]
        # Update the row in the sample table.
        cmd = "UPDATE sample SET user_specimen_id=%s, field_call=%s, depth=%s, percent_missing_data_coral=%s, "
        cmd += "percent_heterozygous_coral=%s, percent_acerv_coral=%s, percent_apalm_coral=%s, bcoral_genet_id=%s, "
        cmd += "registry_id=%s, dna_extraction_method=%s, dna_concentration=%s, colony_location=%s, collection_date=%s "
        cmd += "WHERE affy_id=%s;"
        args = [user_specimen_id, field_call, depth, percent_missing_data_coral, percent_heterozygous_coral,
                percent_acerv_coral, percent_apalm_coral, bcoral_genet_id, registry_id, dna_extraction_method,
                dna_concentration, colony_location, collection_date, affy_id]
        self.update(cmd, args)
        self.flush()
        self.log("Updated the sample table row with affy_id '%s'\n" % affy_id)

    def valid_date(self, val):
        # Date strings must be formated as yyyy-mm-dd.
        items = val.split("-")
        if len(items) != 3:
            return False
        try:
            int(items[0])
            int(items[1])
            int(items[2])
        except Exception:
            return False
        return True


if __name__ == '__main__':
    sdu = StagDatabaseUpdater()
    sdu.run()
    sdu.shutdown()

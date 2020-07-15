#!/usr/bin/env python


import argparse
import datetime
import dateutil.parser
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


def get_value_from_config(config_defaults, value):
    return config_defaults.get(value, None)


def run_command(cmd):
    fstderr, fherr, fstdout, fhout = get_response_buffers()
    proc = subprocess.Popen(args=cmd, stderr=fherr, stdout=fhout, shell=True)
    rc = proc.wait()
    # Check results.
    fherr.close()
    fhout.close()
    check_execution_errors(rc, fstderr, fstdout)


class BloodDatabaseUpdater(object):
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
        self.blood_cell_ids = []

    def connect_db(self):
        url = make_url(self.args.database_connection_string)
        self.log('Attempting to connect to the database...')
        args = url.translate_connect_args(username='user')
        args.update(url.query)
        assert url.get_dialect().name == 'postgresql', 'This script can only be used with PostgreSQL.'
        self.conn = psycopg2.connect(**args)
        self.log("Successfully connected to the database...")

    def convert_date_string_for_database(self, date_string):
        # The value of date_string is %y/%m/%d with
        # the year being 2 digits (yikes!).
        fixed_century = "20%s" % date_string
        fixed_date = fixed_century.replace("/", "-")
        # Convert the string to a format required for
        # inserting into the database.
        database_format = dateutil.parser.parse(fixed_date)
        return str(database_format)

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

    def log(self, msg):
        self.outfh.write("%s\n" % msg)

    def update_blood_cell_table(self, file_path):
        self.log("Updating the blood_cell table...")
        # Columns in the file are:
        # x y z scale rsig
        blood_cell_table_inserts = 0
        with open(file_path, "r") as fh:
            for line in fh:
                line = line.rstrip("\n")
                items = line.split("\t")
                name = items[0]
                x = items[1]
                y = items[2]
                z = items[3]
                scale = items[4]
                rsig = items[5]
                cmd = "INSERT INTO blood_cell VALUES (nextval('blood_cell_id_seq'), %s, %s, %s, %s, %s, %s, %s, %s) RETURNING id;"
                args = ['NOW()', 'NOW()', name, x, y, z, scale, rsig]
                cur = self.update(cmd, args)
                self.flush()
                blood_cell_id = cur.fetchone()[0]
                blood_cell_table_inserts += 1
            self.blood_cell_ids.append(blood_cell_id)
        self.log("Inserted %d rows into the blood table..." % blood_cell_table_inserts)

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
            self.update_blood_cell_table(os.path.join(input_dir, file_name))

    def shutdown(self):
        self.log("Shutting down...")
        self.outfh.flush()
        self.outfh.close()
        self.conn.close()

    def update(self, cmd, args):
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
        return cur


if __name__ == '__main__':
    bdu = BloodDatabaseUpdater()
    bdu.run()
    bdu.shutdown()

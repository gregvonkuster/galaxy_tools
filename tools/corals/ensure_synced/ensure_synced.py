#!/usr/bin/env python

import argparse
import sys

import psycopg2
from sqlalchemy import create_engine
from sqlalchemy import MetaData
from sqlalchemy.engine.url import make_url

metadata = MetaData()

SKIP_VALS = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']


class EnsureSynced(object):
    def __init__(self):
        self.args = None
        self.conn = None
        self.parse_args()
        self.outfh = open(self.args.output, "w")
        self.connect_db()
        self.engine = create_engine(self.args.database_connection_string)
        self.metadata = MetaData(self.engine)
        self.coral_mlg_rep_sample_ids_from_db = []
        self.affy_ids_from_file = []

    def connect_db(self):
        url = make_url(self.args.database_connection_string)
        args = url.translate_connect_args(username='user')
        args.update(url.query)
        assert url.get_dialect().name == 'postgresql', 'This script can only be used with PostgreSQL.'
        self.conn = psycopg2.connect(**args)

    def get_coral_mlg_rep_sample_ids_from_db(self):
        cmd = "SELECT coral_mlg_rep_sample_id, coral_mlg_clonal_id FROM genotype WHERE coral_mlg_rep_sample_id IS NOT NULL AND coral_mlg_rep_sample_id != '' AND coral_mlg_clonal_id != 'failed' ORDER BY coral_mlg_rep_sample_id;"
        cur = self.conn.cursor()
        cur.execute(cmd)
        rows = cur.fetchall()
        for row in rows:
            self.coral_mlg_rep_sample_ids_from_db.append(row[0])
        self.coral_mlg_rep_sample_ids_from_db.sort()

    def get_affy_ids_from_file(self, f):
        with open(f) as fh:
            for line in fh:
                line = line.strip()
                if line in SKIP_VALS:
                    # Skip the first 9 lines in the file.
                    continue
                self.affy_ids_from_file.append(line)
        self.affy_ids_from_file.sort()

    def get_difference(self, list1, list2):
        if len(list1) > len(list2):
            return list(set(list1) - set(list2))
        return list(set(list2) - set(list1))

    def log(self, msg):
        self.outfh.write("%s\n" % msg)

    def parse_args(self):
        parser = argparse.ArgumentParser()
        parser.add_argument('--database_connection_string', dest='database_connection_string', help='Postgres database connection string'),
        parser.add_argument('--affy_ids_from_file', dest='affy_ids_from_file', help='Affy ids taken from all previously genotyped samples vcf file')
        parser.add_argument('--output', dest='output', help='Output dataset'),
        self.args = parser.parse_args()

    def run(self):
        self.get_coral_mlg_rep_sample_ids_from_db()
        self.get_affy_ids_from_file(self.args.affy_ids_from_file)
        if self.coral_mlg_rep_sample_ids_from_db == self.affy_ids_from_file:
            in_sync = True
            self.log("The selected file is in sync with the database.\n\n")
        else:
            in_sync = False
            self.log("The selected file is not in sync with the database.\n\n")
        num_coral_mlg_rep_sample_ids_from_db = len(self.coral_mlg_rep_sample_ids_from_db)
        self.log("Number of coral mlg rep sample ids in the database: %d\n" % num_coral_mlg_rep_sample_ids_from_db)
        num_affy_ids_from_file = len(self.affy_ids_from_file)
        self.log("Number of Affymetrix ids in the file: %d\n" % num_affy_ids_from_file)
        if not in_sync:
            if num_coral_mlg_rep_sample_ids_from_db > num_affy_ids_from_file:
                self.log("The database contains the following Affymetrix ids that are not in the file.\n")
            else:
                self.log("The file contains the following Affymetrix ids that are not in the database.\n")
            diff_list = self.get_difference(self.coral_mlg_rep_sample_ids_from_db, self.affy_ids_from_file)
            for affy_id in diff_list:
                self.log("%s\n" % affy_id)
            self.outfh.flush()
            self.outfh.close()
            sys.exit(1)

    def shutdown(self):
        self.outfh.flush()
        self.outfh.close()
        self.conn.close()


if __name__ == '__main__':
    es = EnsureSynced()
    es.run()
    es.shutdown()

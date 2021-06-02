#!/usr/bin/env python
import argparse
import sys

import psycopg2

from sqlalchemy import MetaData, create_engine
from sqlalchemy.engine.url import make_url


class UniqueMGLIDGenerator(object):
    def __init__(self):
        self.args = None
        self.conn = None
        self.parse_args()
        self.outfh = open(self.args.output, "w")
        self.connect_db()
        self.engine = create_engine(self.args.database_connection_string)
        self.metadata = MetaData(self.engine)

    def parse_args(self):
        parser = argparse.ArgumentParser()
        parser.add_argument('--database_connection_string', dest='database_connection_string', help='Postgres database connection string'),
        parser.add_argument('--output', dest='output', help='Output dataset'),
        self.args = parser.parse_args()

    def connect_db(self):
        url = make_url(self.args.database_connection_string)
        args = url.translate_connect_args(username='user')
        args.update(url.query)
        assert url.get_dialect().name == 'postgresql', 'This script can only be used with PostgreSQL.'
        self.conn = psycopg2.connect(**args)

    def run(self):
        cmd = """
              SELECT DISTINCT coral_mlg_rep_sample_id
              FROM genotype
              WHERE coral_mlg_rep_sample_id is not NULL
              ORDER BY coral_mlg_rep_sample_id;
        """
        cur = self.conn.cursor()
        cur.execute(cmd)
        rows = cur.fetchall()
        for tup in rows:
            self.outfh.write("%s\n" % tup[0])
        self.outfh.close()

    def shutdown(self):
        self.conn.close()

    def stop_err(self, msg):
        sys.stderr.write(msg)
        self.outfh.flush()
        self.outfh.close()
        sys.exit(1)


if __name__ == '__main__':
    umlgidg = UniqueMGLIDGenerator()
    umlgidg.run()
    umlgidg.shutdown()

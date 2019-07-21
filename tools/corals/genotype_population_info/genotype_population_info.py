#!/usr/bin/env python
import argparse
import sys

import psycopg2
from sqlalchemy import create_engine, MetaData
from sqlalchemy.engine.url import make_url


class GenotypeInfoGenerator(object):
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
        parser.add_argument('--input_partial_info', dest='input_partial_info', help='Tabular file containing part of the genotype info')
        parser.add_argument('--output', dest='output', help='Output dataset'),
        self.args = parser.parse_args()

    def connect_db(self):
        url = make_url(self.args.database_connection_string)
        args = url.translate_connect_args(username='user')
        args.update(url.query)
        assert url.get_dialect().name == 'postgresql', 'This script can only be used with PostgreSQL.'
        self.conn = psycopg2.connect(**args)

    def run(self):
        sql = """
              SELECT sample.user_specimen_id,
                     reef.region
              FROM sample
              LEFT OUTER JOIN colony
                              ON sample.colony_id = colony.id
              LEFT OUTER JOIN reef
                              ON reef.id = colony.reef_id
              WHERE sample.affy_id = '%s';
        """
        with open(self.args.input_partial_info, "r") as fh:
            for line in fh:
                line = line.strip()
                out_items = []
                items = line.split('\t')
                # Item number.
                out_items.append(items[0])
                affy_id = items[1]
                out_items.append(affy_id)
                if len(items) == 2:
                    # Example line:
                    # 1 a100000-4368120-060520-256_I07.CEL
                    # The line is missing the user_specimen_id and
                    # region, so retrieve it from the database.
                    query = sql % affy_id
                    cur = self.conn.cursor()
                    cur.execute(query)
                    try:
                        missing_items = cur.fetchone()
                        # user_specimen_id
                        out_items.append(missing_items[0])
                        # region
                        out_items.append(missing_items[1])
                    except Exception as e:
                        msg = "Error retrieving user_specimen_id and region from the database for affy_id %s: %s" % (affy_id, e)
                        self.stop_err(msg)
                else:
                    # The line contains all of the information we need.
                    # user_specimen_id
                    out_items.append(items[3])
                    # region
                    out_items.append(items[9])
                self.outfh.write("%s\n" % "\t".join(out_items))
        self.outfh.close()

    def shutdown(self):
        self.conn.close()

    def stop_err(self, msg):
        sys.stderr.write(msg)
        self.outfh.flush()
        self.outfh.close()
        sys.exit(1)


if __name__ == '__main__':
    gig = GenotypeInfoGenerator()
    gig.run()
    gig.shutdown()

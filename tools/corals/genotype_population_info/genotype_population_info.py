#!/usr/bin/env python
"""
Generate the genotype_population_info.txt file by parsing the information from a VCF
file and querying the stag database that is required to be available within the Galaxy
instance in which this tool is executing.  PostgreSQL 9.1 or greater is required.
"""
import argparse
import sys

import psycopg2
from sqlalchemy import create_engine
from sqlalchemy.engine.url import make_url


class PopInfoGenerator(object):

    def __init__(self):
        self.args = None
        self.conn = None
        self.parse_args()
        self.connect_db()
        self.engine = create_engine(self.args.database_connection_string)

    def parse_args(self):
        parser = argparse.ArgumentParser()
        parser.add_argument('--database_connection_string', dest='database_connection_string', help='Postgres database connection string'),
        parser.add_argument('--input_vcf', dest='input_vcf', help='Input VCF file')
        parser.add_argument('--output', dest='output', help='Output dataset'),
        self.args = parser.parse_args()

    def connect_db(self):
        url = make_url(self.args.database_connection_string)
        self.log('Connecting to database with URL: %s' % url)
        args = url.translate_connect_args(username='user')
        args.update(url.query)
        assert url.get_dialect().name == 'postgresql', 'This script can only be used with PostgreSQL.'
        self.conn = psycopg2.connect(**args)

    def run(self):
        self.gen_pop_info()
        self.fh.flush()
        self.fh.close()

    def shutdown(self):
        self.conn.close()

    def stop_err(self, msg):
        sys.stderr.write(msg)

    def log(self, msg):
        self.fh.write("%s\n" % msg)
        self.fh.flush()

    def get_sample_list(self):
        # Parse the input_vcf file, looking for the first line
        # that starts with the string "#CHROM"
        with open(self.args.input_vcf, "r") as vcfh:
            for line in vcfh:
                if not line.startswith("#CHROM"):
                    continue
                line = line.rstrip("\r\n")
                # Example line:
                # #CHROM  13704   13706   13708   13736   13748   13762   13782
                items = line.split("\t")
                sample_list = items[8:]
                break
        return sample_list

    def get_region_list(self, sample_list):
        # Retrieve the value of the region column in the reef table
        # for each sample_id in the sample_list.
        region_list = []
        for sample_id in sample_list:
            sql = """SELECT reef.region
                  FROM reef
                  LEFT OUTER JOIN colony ON reef.id = colony.reef_id
                  LEFT OUTER JOIN sample ON sample.colony_id = colony.id
                  WHERE sample.id = '%s';""" % sample_id
            cur = self.conn.cursor()
            cur.execute(sql)
            region_list.append(cur.fetchone()[0])
        return region_list

    def gen_pop_info(self):
        sample_list = self.get_sample_list()
        region_list = self.get_region_list(sample_list)
        # The output file will consist of columns:
        # Item # Sample ID Region
        with open(self.args.output, "w") as outfh:
            for i, sample_id in sample_list:
                outfh.write("%d\t%s\t%s\n" % (i, sample_id, region_list[1]))


if __name__ == '__main__':
    pop_info_generator = PopInfoGenerator()
    pop_info_generator.run()
    pop_info_generator.shutdown()

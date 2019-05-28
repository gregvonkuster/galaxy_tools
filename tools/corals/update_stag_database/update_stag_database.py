from __future__ import print_function

import argparse
import datetime
import dateutil.parser
import psycopg2
import sys

from sqlalchemy import create_engine
from sqlalchemy import MetaData
from sqlalchemy import sql
from sqlalchemy.engine.url import make_url

now = datetime.datetime.utcnow
metadata = MetaData()


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


def handle_null(column_val):
    if set_to_null(column_val):
        return sql.null(), "%s"
    else:
        return column_val, "'%s'"


def set_to_null(val):
    if val in ["", "NA", "NULL"]:
        return True
    return False


class StagDatabaseUpdater(object):
    def __init__(self):
        self.args = None
        self.conn = None
        self.parse_args()
        self.year_from_now = get_year_from_now()
        self.outfh = open(self.args.output, "w")
        self.connect_db()
        self.engine = create_engine(self.args.database_connection_string)
        self.metadata = MetaData(self.engine)
        self.genotype_ids = []
        self.genotype_table_inserts = 0
        self.person_ids = []
        self.person_table_inserts = 0
        self.phenotype_ids = []
        self.phenotype_table_inserts = 0
        self.reef_ids = []
        self.reef_table_inserts = 0
        self.sample_table_inserts = 0
        self.taxonomy_ids = []
        self.taxonomy_table_inserts = 0

    def connect_db(self):
        url = make_url(self.args.database_connection_string)
        self.log('Connecting to database with URL: %s' % url)
        args = url.translate_connect_args(username='user')
        args.update(url.query)
        assert url.get_dialect().name == 'postgresql', 'This script can only be used with PostgreSQL.'
        self.conn = psycopg2.connect(**args)

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

    def load_genotype_table(self, file_path):
        # Columns in the genotype file are:
        # affy_id coral_mlg_clonal_id.x coral_mlg_rep_sample_id.x symbio_mlg_clonal_id symbio_mlg_rep_sample_id
        # genetic_coral_species_call bcoral_genet_id bsym_genet_id DB_match coral_mlg_clonal_id.y
        # coral_mlg_rep_sample_id.y
        with open(file_path) as fh:
            for i, line in enumerate(fh):
                if i == 0:
                    # Skip header
                    continue
                items = line.split("\t")
                coral_mlg_clonal_id = items[1]
                # In those cases in which column 2 is null, we'll set
                # coral_mlg_rep_sample_id to the value of column 10.
                if set_to_null(items[2]):
                    coral_mlg_rep_sample_id = items[10]
                else:
                    coral_mlg_rep_sample_id = items[2]
                coral_mlg_rep_sample_id, coral_mlg_rep_sample_id_db_str = handle_null(coral_mlg_rep_sample_id)
                symbio_mlg_clonal_id, symbio_mlg_clonal_id_db_str = handle_null(items[3])
                symbio_mlg_rep_sample_id, symbio_mlg_rep_sample_id_db_str = handle_null(items[4])
                genetic_coral_species_call, genetic_coral_species_call_db_str = handle_null(items[5])
                bcoral_genet_id, bcoral_genet_id_db_str = handle_null(items[6])
                bsym_genet_id, bsym_genet_id_db_str = handle_null(items[7])
                # See if we need to add a row to the table.
                # TODO: find out if the coral_mlg_clonal_id column is the
                # optimal unique identifier for determining if a new row
                # should be inserted.
                cmd = "SELECT id FROM genotype WHERE coral_mlg_clonal_id = '%s';" % coral_mlg_clonal_id
                cur = self.conn.cursor()
                cur.execute(cmd)
                try:
                    genotype_id = cur.fetchone()[0]
                except Exception:
                    # Insert a row into the genotype table.
                    cmd = "INSERT INTO genotype VALUES (nextval('genotype_id_seq'), '%s' %s %s %s %s %s %s) RETURNING id;"
                    cmd = cmd % (coral_mlg_clonal_id,
                                 coral_mlg_rep_sample_id_db_str,
                                 symbio_mlg_clonal_id_db_str,
                                 symbio_mlg_rep_sample_id_db_str,
                                 genetic_coral_species_call_db_str,
                                 bcoral_genet_id_db_str,
                                 bsym_genet_id_db_str)
                    args = [coral_mlg_rep_sample_id, symbio_mlg_clonal_id, symbio_mlg_rep_sample_id, genetic_coral_species_call, bcoral_genet_id, bsym_genet_id]
                    cur = self.update(cmd, tuple(args))
                    self.flush()
                    genotype_id = cur.fetchone()[0]
                    self.genotype_table_inserts += 1
                self.genotype_ids.append(genotype_id)

    def load_person_table(self, file_path):
        # Columns in the person file are:
        # last_name first_name organization email
        with open(file_path) as fh:
            for i, line in enumerate(fh):
                if i == 0:
                    # Skip header
                    continue
                items = line.split("\t")
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
                    cmd = "INSERT INTO person VALUES (nextval('person_id_seq'), '%s', '%s', '%s', '%s') RETURNING id;"
                    args = [last_name, first_name, organization, email]
                    cur = self.update(cmd, tuple(args))
                    self.flush()
                    person_id = cur.fetchone()[0]
                    self.person_table_inserts += 1
                self.person_ids.append(person_id)

    def load_phenotype_table(self, file_path):
        # Columns in the phenotype file are:
        # last_name first_name organization email disease_resist
        # bleach_resist mortality tle spawning sperm_motility
        # healing_time
        with open(file_path) as fh:
            for i, line in enumerate(fh):
                if i == 0:
                    # Skip header
                    continue
                items = line.split("\t")
                disease_resist, disease_resist_db_str = handle_null(items[4])
                bleach_resist, bleach_resist_db_str = handle_null(items[5])
                mortality, mortality_db_str = handle_null(items[6])
                tle, tle_db_str = handle_null(items[7])
                spawning, spawning_db_str = handle_null(items[8])
                sperm_motility, sperm_motility_db_str = handle_null(items[9])
                healing_time, healing_time_db_str = handle_null(items[10])
                # See if we need to add a row to the phenotype table.
                cmd = """
                    SELECT id FROM phenotype
                    WHERE disease_resist = %s
                    AND bleach_resist = %s
                    AND mortality = %s
                    AND tle = %s
                    AND spawning = %s
                    AND sperm_motility = %s
                    AND healing_time = %s;""" % (disease_resist_db_str,
                                                 bleach_resist_db_str,
                                                 mortality_db_str,
                                                 tle_db_str,
                                                 spawning_db_str,
                                                 sperm_motility_db_str,
                                                 healing_time_db_str)
                cur = self.conn.cursor()
                cur.execute(cmd)
                try:
                    phenotype_id = cur.fetchone()[0]
                except Exception:
                    # Insert a row into the phenotype table.
                    cmd = "INSERT INTO phenotype VALUES (nextval('phenotype_id_seq'), %s %s %s %s %s %s %s) RETURNING id;" % \
                        (disease_resist_db_str, bleach_resist_db_str, mortality_db_str, tle_db_str, spawning_db_str, sperm_motility_db_str, healing_time_db_str)
                    args = [disease_resist, bleach_resist, mortality, tle, spawning, sperm_motility, healing_time]
                    cur = self.update(cmd, tuple(args))
                    self.flush()
                    phenotype_id = cur.fetchone()[0]
                    self.phenotype_table_inserts += 1
                self.phenotype_ids.append(phenotype_id)

    def load_reef_table(self, file_path):
        # Columns in the reef file are:
        # name region latitude longitude geographic_origin
        with open(file_path) as fh:
            for i, line in enumerate(fh):
                if i == 0:
                    # Skip header
                    continue
                items = line.split("\t")
                name = items[0]
                region = items[1]
                latitude = "%6f" % float(items[2])
                longitude = "%6f" % float(items[3])
                geographic_origin = items[5]
                # See if we need to add a row to the reef table.
                cmd = "SELECT id FROM reef WHERE name = '%s' AND region = '%s' AND latitude = %s AND longitude = %s AND geographic_origin = '%s';" % \
                    (name, region, latitude, longitude, geographic_origin)
                cur = self.conn.cursor()
                cur.execute(cmd)
                try:
                    reef_id = cur.fetchone()[0]
                except Exception:
                    # Insert a row into the reef table.
                    cmd = "INSERT INTO reef VALUES (nextval('reef_id_seq'), '%s' '%s' %s %s '%s') RETURNING id;"
                    args = [name, region, latitude, longitude, geographic_origin]
                    cur = self.update(cmd, tuple(args))
                    self.flush()
                    reef_id = cur.fetchone()[0]
                    self.reef_table_inserts += 1
                self.reef_ids.append(reef_id)

    def load_sample_table(self, file_path):
        # Columns in the sample file are:
        # affy_id colony_location collection_date user_specimen_id registry_id
        # depth dna_extraction_method dna_concentration public public_after_date
        # percent_missing_data_coral percent_missing_data_sym percent_reference_coral percent_reference_sym percent_alternative_coral
        # percent_alternative_sym percent_heterozygous_coral percent_heterozygous_sym field_call
        with open(file_path) as fh:
            for i, line in enumerate(fh):
                if i == 0:
                    # Skip header
                    continue
                # Keep track of foreign keys since we skip the header line
                id_index = i - 1
                items = line.split("\t")
                affy_id = items[0]
                sample_id = self.get_next_sample_id()
                # FIXME: need to insert alleles soe we have the list of allele_ids.
                allele_id = sql.null()
                genotype_id = self.genotype_ids[id_index]
                phenotype_id = self.phenotype_ids[id_index]
                # FIXME: We cannot populate the experiment table with our current data.
                experiment_id = sql.null()

    def load_taxonomy_table(self, file_path):
        # Columns in the taxonomy file are:
        # genetic_coral_species_call affy_id genus_name species_name"
        with open(file_path) as fh:
            for i, line in enumerate(fh):
                if i == 0:
                    # Skip header
                    continue
                items = line.split("\t")
                genus_name = items[2]
                species_name = items[3]
                # See if we need to add a row to the taxonomy table.
                cmd = "SELECT id FROM taxonomy WHERE species_name = '%s' AND genus_name = '%s';" % (species_name, genus_name)
                cur = self.conn.cursor()
                cur.execute(cmd)
                try:
                    taxonomy_id = cur.fetchone()[0]
                except Exception:
                    # Insert a row into the taxonomy table.
                    cmd = "INSERT INTO reef VALUES (nextval('taxonomy_id_seq'), '%s' '%s') RETURNING id;"
                    args = [species_name, genus_name]
                    cur = self.update(cmd, tuple(args))
                    self.flush()
                    taxonomy_id = cur.fetchone()[0]
                    self.taxonomy_table_inserts += 1
                self.taxonomy_ids.append(taxonomy_id)

    def parse_args(self):
        parser = argparse.ArgumentParser()
        parser.add_argument('--database_connection_string', dest='database_connection_string', help='Postgres database connection string'),
        parser.add_argument('--input', dest='inputs', action='append', nargs=6, help='Input datasets for database insertion')
        parser.add_argument('--output', dest='output', help='Output dataset'),
        self.args = parser.parse_args()

    def run(self):
        for input in self.args.inputs:
            # Tables must be loaded in such a way that foreign keys
            # are properly handled.  The sample table must be loaded
            # last.
            file_path, file_name = input
            if file_name.startswith("genotype"):
                genotype_file = file_path
            elif file_name.startswith("person"):
                person_file = file_path
            elif file_name.startswith("phenotype"):
                phenotype_file = file_path
            elif file_name.startswith("reef"):
                reef_file = file_path
            elif file_name.startswith("sample"):
                sample_file = file_path
            elif file_name.startswith("taxonomy"):
                taxonomy_file = file_path
        # Now tables can be loaded in the appropriate order.
        self.load_genotype_table(genotype_file)
        self.load_person_table(person_file)
        self.load_phenotype_table(phenotype_file)
        self.load_reef_table(reef_file)
        self.load_taxonomy_table(taxonomy_file)
        self.load_sample_table(sample_file)

    def shutdown(self):
        self.conn.close()

    def stop_err(self, msg):
        sys.stderr.write(msg)
        self.outfh.flush()
        self.outfh.close()
        sys.exit(1)

    def update(self, sql, args):
        try:
            cur = self.conn.cursor()
            cur.execute(sql, args)
        except Exception as e:
            msg = "Caught exception executing SQL:\n%s\nException:\n%s\n" % (sql.format(args), e)
            self.stop_err(msg)
        return cur


if __name__ == '__main__':
    sdu = StagDatabaseUpdater()
    sdu.run()
    sdu.shutdown()

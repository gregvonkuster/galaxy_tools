from __future__ import print_function

import argparse
import datetime
import dateutil.parser
import os
import psycopg2
import sys

from sqlalchemy import create_engine
from sqlalchemy import MetaData
from sqlalchemy import sql
from sqlalchemy.engine.url import make_url

now = datetime.datetime.utcnow
metadata = MetaData()


def get_sql_param_val_str(column_val):
    if set_to_null(column_val):
        return "is null"
    else:
        return "= '%s'" % column_val


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


def handle_null(val):
    if set_to_null(val):
        return None
    return val


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
        self.colony_ids = []
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

    def update_colony_table(self, file_path):
        self.log("Updating the colony table...")
        # Columns in the colony file are:
        # latitude longitude depth
        colony_table_inserts = 0
        with open(file_path) as fh:
            for i, line in enumerate(fh):
                if i == 0:
                    # Skip header
                    continue
                line = line.rstrip()
                items = line.split("\t")
                latitude = items[0]
                latitude_param_val_str = get_sql_param_val_str(latitude)
                longitude = items[1]
                longitude_param_val_str = get_sql_param_val_str(longitude)
                depth = items[2]
                depth_param_val_str = get_sql_param_val_str(depth)
                # See if we need to add a row to the table.
                cmd = "SELECT id FROM colony WHERE latitude %s " % latitude_param_val_str
                cmd += "AND longitude %s AND depth %s;" % (longitude_param_val_str, depth_param_val_str)
                cur = self.conn.cursor()
                cur.execute(cmd)
                try:
                    colony_id = cur.fetchone()[0]
                except Exception:
                    # Insert a row into the colony table.
                    cmd = "INSERT INTO colony VALUES (nextval('colony_id_seq'), %s, %s, "
                    cmd += "$$%s$$, $$%s$$, $$%s$$) RETURNING id;"
                    args = ['NOW()', 'NOW()', longitude, latitude, depth]
                    cur = self.update(cmd, args)
                    self.flush()
                    colony_id = cur.fetchone()[0]
                    colony_table_inserts += 1
                self.colony_ids.append(colony_id)
        self.log("Inserted %d rows into the colony table..." % colony_table_inserts)

    def update_genotype_table(self, file_path):
        self.log("Updating the genotype table...")
        # Columns in the genotype file are:
        # affy_id coral_mlg_clonal_id.x coral_mlg_rep_sample_id.x symbio_mlg_clonal_id symbio_mlg_rep_sample_id
        # genetic_coral_species_call bcoral_genet_id bsym_genet_id DB_match coral_mlg_clonal_id.y
        # coral_mlg_rep_sample_id.y
        genotype_table_inserts = 0
        with open(file_path) as fh:
            for i, line in enumerate(fh):
                if i == 0:
                    # Skip header
                    continue
                line = line.rstrip()
                items = line.split("\t")
                coral_mlg_clonal_id = items[1]
                # In those cases in which column 2 is null, we'll set
                # coral_mlg_rep_sample_id to the value of column 10.
                if set_to_null(items[2]):
                    coral_mlg_rep_sample_id = items[10]
                else:
                    coral_mlg_rep_sample_id = items[2]
                coral_mlg_rep_sample_id_param_val_str = get_sql_param_val_str(coral_mlg_rep_sample_id)
                symbio_mlg_clonal_id = items[3]
                symbio_mlg_clonal_id_param_val_str = get_sql_param_val_str(symbio_mlg_clonal_id)
                symbio_mlg_rep_sample_id = items[4]
                symbio_mlg_rep_sample_id_param_val_str = get_sql_param_val_str(symbio_mlg_rep_sample_id)
                genetic_coral_species_call = items[5]
                genetic_coral_species_call_param_val_str = get_sql_param_val_str(genetic_coral_species_call)
                bcoral_genet_id = items[6]
                bcoral_genet_id_param_val_str = get_sql_param_val_str(bcoral_genet_id)
                bsym_genet_id = items[7]
                bsym_genet_id_param_val_str = get_sql_param_val_str(bsym_genet_id)
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
                    cmd = "INSERT INTO genotype VALUES (nextval('genotype_id_seq'), NOW(), NOW(), "
                    cmd += "$$%s$$, $$%s$$, $$%s$$, $$%s$$, $$%s$$, $$%s$$, $$%s$$) RETURNING id;"
                    cmd = cmd % (coral_mlg_clonal_id, coral_mlg_rep_sample_id, symbio_mlg_clonal_id,
                                 symbio_mlg_rep_sample_id, genetic_coral_species_call, bcoral_genet_id,
                                 bsym_genet_id)
                    args = [coral_mlg_rep_sample_id, symbio_mlg_clonal_id, symbio_mlg_rep_sample_id,
                            genetic_coral_species_call, bcoral_genet_id, bsym_genet_id]
                    cur = self.update(cmd, args)
                    self.flush()
                    genotype_id = cur.fetchone()[0]
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
                line = line.strip()
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
                    cmd = "INSERT INTO person VALUES (nextval('person_id_seq'), NOW(), NOW(), "
                    cmd += "$$%s$$, $$%s$$, $$%s$$, $$%s$$) RETURNING id;"
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
        # last_name first_name organization email disease_resist
        # bleach_resist mortality tle spawning sperm_motility
        # healing_time
        phenotype_table_inserts = 0
        with open(file_path) as fh:
            for i, line in enumerate(fh):
                if i == 0:
                    # Skip header
                    continue
                line = line.rstrip()
                items = line.split("\t")
                disease_resist = items[4]
                disease_resist_param_val_str = get_sql_param_val_str(disease_resist)
                bleach_resist = items[5]
                bleach_resist_param_val_str = get_sql_param_val_str(bleach_resist)
                mortality = items[6]
                mortality_param_val_str = get_sql_param_val_str(mortality)
                tle = items[7]
                tle_param_val_str = get_sql_param_val_str(tle)
                spawning = items[8]
                spawning_param_val_str = get_sql_param_val_str(spawning)
                sperm_motility = items[9]
                sperm_motility_param_val_str = get_sql_param_val_str(sperm_motility)
                healing_time = items[10]
                healing_time_param_val_str = get_sql_param_val_str(healing_time)
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
                    cmd += "$$%s$$, $$%s$$, $$%s$$, " % (disease_resist, bleach_resist, mortality)
                    cmd += "$$%s$$, $$%s$$, $$%s$$, " % (tle, spawning, sperm_motility)
                    cmd += "$$%s$$) RETURNING id;" % healing_time
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
        reef_table_inserts = 0
        with open(file_path) as fh:
            for i, line in enumerate(fh):
                if i == 0:
                    # Skip header
                    continue
                line = line.rstrip()
                items = line.split("\t")
                name = items[0]
                region = items[1]
                latitude = "%6f" % float(items[2])
                longitude = "%6f" % float(items[3])
                geographic_origin = items[4]
                # See if we need to add a row to the reef table.
                cmd = "SELECT id FROM reef WHERE name = '%s' AND region = '%s' " % (name, region)
                cmd += "AND latitude = %s AND longitude = %s " % (latitude, longitude)
                cmd += "AND geographic_origin = '%s';" % geographic_origin
                cur = self.conn.cursor()
                cur.execute(cmd)
                try:
                    reef_id = cur.fetchone()[0]
                except Exception:
                    # Insert a row into the reef table.
                    cmd = "INSERT INTO reef VALUES (nextval('reef_id_seq'), %s, %s, "
                    cmd += "$$%s$$, $$%s$$, $$%s$$, $$%s$$, $$%s$$) RETURNING id;"
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
        # percent_alternative_sym percent_heterozygous_coral percent_heterozygous_sym field_call
        sample_table_inserts = 0
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
                # FIXME: need to insert alleles so we have the list of allele_ids.
                allele_id = sql.null()
                genotype_id = self.genotype_ids[id_index]
                phenotype_id = self.phenotype_ids[id_index]
                # FIXME: We cannot populate the experiment table with our current data.
                experiment_id = sql.null()
                colony_id = self.colony_ids[id_index]
                colony_location, colony_location_param_val_str = get_sql_param_val_str(items[1])
                # FIXME: We cannot populate the fragment table with our current data.
                fragment_id = sql.null()
                taxonomy_id = self.taxonomy_ids[id_index]
                # FIXME: We should eliminate the collector table and use the person table.
                collector_id = self.person_ids[id_index]
                collection_date = items[2]
                user_specimen_id = items[3]
                registry_id, registry_id_param_val_str = get_sql_param_val_str(items[4])
                depth = items[5]
                dna_extraction_method, dna_extraction_method_param_val_str = get_sql_param_val_str(items[6])
                dna_concentration = items[7]
                dna_concentration_param_val_str = get_sql_param_val_str(dna_concentration)
                public = items[8]
                if set_to_null(items[9]):
                    public_after_date = self.year_from_now
                else:
                    public_after_date = items[9]
                percent_missing_data_coral = items[10]
                percent_missing_data_sym = items[11]
                percent_reference_coral = items[12]
                percent_reference_sym = items[13]
                percent_alternative_coral = items[14]
                percent_alternative_sym = items[15]
                percent_heterozygous_coral = items[16]
                percent_heterozygous_sym = items[17]
                field_call = items[18]
                # Insert a row into the sample table.
                cmd = """
                    INSERT INTO sample
                    VALUES (nextval('sample_id_seq'),
                            '%s' '%s' %s %s %s, %s, %s, '%s', %s, %s, %s, '%s', '%s', '%s', %s, '%s', '%s', '%s', '%s', %s, %s, %s, %s, %s, %s, %s, %s, '%s')
                    RETURNING id;"""
                args = [affy_id, sample_id, allele_id, genotype_id, phenotype_id,
                        experiment_id, colony_id, colony_location, fragment_id, taxonomy_id,
                        collector_id, collection_date, user_specimen_id, registry_id, depth,
                        dna_extraction_method, dna_concentration, public, public_after_date,
                        percent_missing_data_coral, percent_missing_data_sym, percent_reference_coral,
                        percent_reference_sym, percent_alternative_coral, percent_alternative_sym,
                        percent_heterozygous_coral, percent_heterozygous_sym, field_call]
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
                    cmd = "INSERT INTO taxonomy VALUES (nextval('taxonomy_id_seq'), '%s' '%s') RETURNING id;"
                    args = [species_name, genus_name]
                    cur = self.update(cmd, args)
                    self.flush()
                    taxonomy_id = cur.fetchone()[0]
                    taxonomy_table_inserts += 1
                self.taxonomy_ids.append(taxonomy_id)
        self.log("Inserted %d rows into the taxonomy table..." % taxonomy_table_inserts)

    def parse_args(self):
        parser = argparse.ArgumentParser()
        parser.add_argument('--database_connection_string', dest='database_connection_string', help='Postgres database connection string'),
        parser.add_argument('--input_dir', dest='input_dir', help='Input datasets for database insertion')
        parser.add_argument('--output', dest='output', help='Output dataset'),
        self.args = parser.parse_args()

    def run(self):
        input_dir = self.args.input_dir
        for file_name in os.listdir(input_dir):
            # Tables must be loaded in such a way that foreign keys
            # are properly handled.  The sample table must be loaded
            # last.
            if file_name.startswith("colony"):
                colony_file = os.path.join(input_dir, file_name)
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
        self.update_genotype_table(genotype_file)
        self.update_person_table(person_file)
        self.update_phenotype_table(phenotype_file)
        self.update_reef_table(reef_file)
        self.update_colony_table(colony_file)
        #self.update_taxonomy_table(taxonomy_file)
        #self.update_sample_table(sample_file)

    def shutdown(self):
        self.log("Shutting down...")
        self.conn.close()

    def stop_err(self, msg):
        sys.stderr.write(msg)
        self.outfh.flush()
        self.outfh.close()
        sys.exit(1)

    def update(self, sql, args):
        for i, arg in enumerate(args):
            args[i] = handle_null(arg)
        try:
            cur = self.conn.cursor()
            cur.execute(sql, tuple(args))
        except Exception as e:
            msg = "Caught exception executing SQL:\n%s\nException:\n%s\n" % (sql.format(args), e)
            self.stop_err(msg)
        return cur


if __name__ == '__main__':
    sdu = StagDatabaseUpdater()
    sdu.run()
    sdu.shutdown()

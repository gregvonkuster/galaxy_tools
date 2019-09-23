#!/usr/bin/env python
from __future__ import print_function

import argparse
import datetime
import psycopg2

from sqlalchemy import create_engine
from sqlalchemy import MetaData
from sqlalchemy.engine.url import make_url

metadata = MetaData()

COLUMNS = ["Affymetrix ID", "Sample ID", "User Specimen ID", "Field Call", "Sample Depth",
           "Percent Missing Data Coral", "Percent Heterozygous Coral", "Percent Acerv Coral",
           "Percent Apalm Coral", "Bcoral Genet Id", "Registry ID", "DNA Extraction Method",
           "DNA Concentration", "Colony Location", "Colony Latitude", "Colony Longitude",
           "Colony Depth", "Reef Name", "Region", "Reef Latitude", "Reef Longitude",
           "GPS Coordinates Associated With", "Coral Mlg Clonal ID", "Coral Mlg Rep Sample ID",
           "Genetic Coral Species Call", "Spawning", "Sperm Motility", "TLE", "Disease Resist",
           "Bleach Resist", "Mortality", "Healing Time", "Sequencing Facility", "Array Version",
           "Plate Barcode", "Collector Last Name", "First Name", "Organization", "Email",
           "Collection Date"]


class ExportAllSampleData(object):
    def __init__(self):
        self.args = None
        self.conn = None
        self.parse_args()
        self.outfh = open(self.args.output, "w")
        self.outfh.write("%s\n" % "\t".join(COLUMNS))
        self.connect_db()
        self.engine = create_engine(self.args.database_connection_string)

    def connect_db(self):
        url = make_url(self.args.database_connection_string)
        args = url.translate_connect_args(username='user')
        args.update(url.query)
        assert url.get_dialect().name == 'postgresql', 'This script can only be used with PostgreSQL.'
        self.conn = psycopg2.connect(**args)

    def export_from_db(self):
        today = datetime.date.today()
        cmd = """
           SELECT sample.affy_id, sample.sample_id, sample.genotype_id,
           sample.phenotype_id, sample.experiment_id, sample.colony_id,
           sample.colony_location, sample.collector_id, sample.collection_date,
           sample.user_specimen_id, sample.registry_id, sample.depth AS sample_depth,
           sample.dna_extraction_method, sample.dna_concentration,
           sample.percent_missing_data_coral, sample.percent_acerv_coral,
           sample.percent_apalm_coral, sample.percent_heterozygous_coral,
           sample.field_call, sample.bcoral_genet_id, genotype.coral_mlg_clonal_id,
           genotype.coral_mlg_rep_sample_id, genotype.genetic_coral_species_call,
           phenotype.spawning, phenotype.sperm_motility, phenotype.tle,
           phenotype.disease_resist, phenotype.bleach_resist, phenotype.mortality,
           phenotype.healing_time, experiment.seq_facility, experiment.array_version,
           experiment.plate_barcode, colony.latitude AS colony_latitude,
           colony.longitude AS colony_longitude, colony.depth AS colony_depth,
           reef.name, reef.region, reef.latitude AS reef_latitude, reef.longitude AS reef_longitude,
           reef.geographic_origin, person.last_name, person.first_name,
           person.organization, person.email
           FROM sample
           LEFT OUTER JOIN genotype
                           ON sample.genotype_id = genotype.id
           LEFT OUTER JOIN phenotype
                           ON sample.phenotype_id = phenotype.id
           LEFT OUTER JOIN experiment
                           ON sample.experiment_id = experiment.id
           LEFT OUTER JOIN colony
                           ON sample.colony_id = colony.id
           LEFT OUTER JOIN reef
                           ON reef.id = colony.reef_id
           LEFT OUTER JOIN person
                           ON sample.collector_id = person.id
           WHERE sample.public OR sample.public_after_date < date'%s'
           ORDER BY affy_id;""" % today
        # Instantiate the cursor.
        cur = self.conn.cursor()
        # Execute the query.
        cur.execute(cmd)
        rows = cur.fetchall()
        for tup in rows:
            values = self.extract_values(tup)
            # Output the row.
            self.outfh.write("%s\n" % "\t".join(values))

    def extract_values(self, tup):
        values = []
        # Extract the items from the tuple.
        affy_id = self.get_value(tup[0])
        sample_id = self.get_value(tup[1])
        colony_location = self.get_value(tup[6])
        collection_date = self.get_value(tup[8])
        if len(collection_date) > 0:
            collection_date = collection_date[:10]
        user_specimen_id = self.get_value(tup[9])
        registry_id = self.get_value(tup[10])
        sample_depth = self.get_value(tup[11])
        dna_extraction_method = self.get_value(tup[12])
        dna_concentration = self.get_value(tup[13])
        percent_missing_data_coral = self.get_value(tup[14])
        percent_acerv_coral = self.get_value(tup[15])
        percent_apalm_coral = self.get_value(tup[16])
        percent_heterozygous_coral = self.get_value(tup[17])
        field_call = self.get_value(tup[18])
        bcoral_genet_id = self.get_value(tup[19])
        coral_mlg_clonal_id = self.get_value(tup[20])
        coral_mlg_rep_sample_id = self.get_value(tup[21])
        genetic_coral_species_call = self.get_value(tup[22])
        spawning = self.get_value(tup[23])
        sperm_motility = self.get_value(tup[24])
        tle = self.get_value(tup[25])
        disease_resist = self.get_value(tup[26])
        bleach_resist = self.get_value(tup[27])
        mortality = self.get_value(tup[28])
        healing_time = self.get_value(tup[29])
        seq_facility = self.get_value(tup[30])
        array_version = self.get_value(tup[31])
        plate_barcode = self.get_value(tup[32])
        colony_latitude = self.get_value(tup[33])
        colony_longitude = self.get_value(tup[34])
        colony_depth = self.get_value(tup[35])
        reef_name = self.get_value(tup[36])
        region = self.get_value(tup[37])
        reef_latitude = self.get_value(tup[38])
        reef_longitude = self.get_value(tup[39])
        geographic_origin = self.get_value(tup[40])
        last_name = self.get_value(tup[41])
        first_name = self.get_value(tup[42])
        organization = self.get_value(tup[43])
        email = self.get_value(tup[44])
        # Append the columns in the specified order.
        values.append(affy_id)
        values.append(sample_id)
        values.append(user_specimen_id)
        values.append(field_call)
        values.append(sample_depth)
        values.append(percent_missing_data_coral)
        values.append(percent_heterozygous_coral)
        values.append(percent_acerv_coral)
        values.append(percent_apalm_coral)
        values.append(bcoral_genet_id)
        values.append(registry_id)
        values.append(dna_extraction_method)
        values.append(dna_concentration)
        values.append(colony_location)
        values.append(colony_latitude)
        values.append(colony_longitude)
        values.append(colony_depth)
        values.append(reef_name)
        values.append(region)
        values.append(reef_latitude)
        values.append(reef_longitude)
        values.append(geographic_origin)
        values.append(coral_mlg_clonal_id)
        values.append(coral_mlg_rep_sample_id)
        values.append(genetic_coral_species_call)
        values.append(spawning)
        values.append(sperm_motility)
        values.append(tle)
        values.append(disease_resist)
        values.append(bleach_resist)
        values.append(mortality)
        values.append(healing_time)
        values.append(seq_facility)
        values.append(array_version)
        values.append(plate_barcode)
        values.append(last_name)
        values.append(first_name)
        values.append(organization)
        values.append(email)
        values.append(collection_date)
        return values

    def get_value(self, loc):
        return str(loc) or ""

    def parse_args(self):
        parser = argparse.ArgumentParser()
        parser.add_argument('--database_connection_string', dest='database_connection_string', help='Postgres database connection string'),
        parser.add_argument('--output', dest='output', help='Output dataset'),
        self.args = parser.parse_args()

    def run(self):
        self.export_from_db()

    def shutdown(self):
        self.outfh.flush()
        self.outfh.close()
        self.conn.close()


if __name__ == '__main__':
    easd = ExportAllSampleData()
    easd.run()
    easd.shutdown()

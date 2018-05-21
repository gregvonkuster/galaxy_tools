#!/usr/bin/env python
"""
add_plant_tribes_scaffold.py - A script for adding a scaffold to the Galaxy PlantTribes
database efficiently by bypassing the Galaxy model and operating directly on the database.
PostgreSQL 9.1 or greater is required.
"""
from __future__ import print_function

import argparse
import glob
import os
import sys

import psycopg2
from sqlalchemy.engine.url import make_url

galaxy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
sys.path.insert(1, os.path.join(galaxy_root, 'lib'))

import galaxy.config
from galaxy.util.script import app_properties_from_args, populate_config_args


class AddScaffold(object):
    def __init__(self):
        self.args = None
        self.config = None
        self.conn = None
        self.scaffold_recs = []
        self.clustering_methods = []
        self.scaffold_genes_dict = {}
        self.species_ids_dict = {}
        self.species_genes_dict = {}
        self.gene_sequences_dict = {}
        self.__parse_args()
        self.__load_config()
        self.__connect_db()

    def __parse_args(self):
        parser = argparse.ArgumentParser()
        populate_config_args(parser)
        parser.add_argument('--dry_run', action='store_true', dest='dry_run', help="Dry run (rollback all transactions)", default=False)
        parser.add_argument('--scaffold_id', dest='scaffold_id', help='Full path to PlantTribes scaffold directory', default=os.path.join(galaxy_root, 'tool-data', 'plant_tribes', 'scaffolds', '22Gv1.1'))
        parser.add_argument('--taxa_lineage', dest='taxa_lineage', help='Scaffold taxa lineage', default=os.path.join(galaxy_root, 'plant_tribes', 'scaffolds', '22Gv1.1', '22Gv1.1taxaLineage'))
        self.args = parser.parse_args()

    def __load_config(self):
        app_properties = app_properties_from_args(self.args)
        self.config = galaxy.config.Configuration(**app_properties)

    def __connect_db(self):
        url = make_url(galaxy.config.get_plant_tribes_database_url(self.config))
        print('Connecting to database with URL: %s' % url)
        args = url.translate_connect_args(username='user')
        args.update(url.query)

        assert url.get_dialect().name == 'postgresql', 'This script can only be used with PostgreSQL.'

        self.conn = psycopg2.connect(**args)

    def _run(self):
        self.process_annot_dir()
        self.process_scaffold_config_files()
        self.process_orthogroup_fasta_files()

    def _update(self, sql, args):
        try:
            cur = self.conn.cursor()
            cur.execute(sql, args)
        except:
            print('Caught exception executing SQL:\n%s\n' % sql.format(args))
            raise
        return cur

    def process_orthogroup_fasta_files(self):
        scaffold_id = os.path.basename(self.args.scaffold_id)
        aa_dict = {}
        dna_dict = {}
        # Populate aa_dict and dna_dict.
        for clustering_method in self.clustering_methods:
            file_dir = os.path.join(self.args.scaffold_id, 'fasta', clustering_method)
            for file_name in os.listdir(file_dir):
                items = file_name.split(".")
                orthogroup_id = items[0]
                file_extension = items[1]
                if file_extension == "fna":
                    adict = dna_dict
                else:
                    adict = aa_dict
                file_path = os.path.join(file_dir, file_name)
                with open(file_path, "r") as fh:
                    for i, line in enumerate(fh):
                        line = line.strip()
                        if len(line) == 0:
                            # Skip blank lines (shoudn't happen).
                            continue
                        if line.startswith(">"):
                            # Example line:
                            # >gnl_Ambtr1.0.27_AmTr_v1.0_scaffold00001.110
                            gene_id = line.lstrip(">")
                            # The dictionary keys will combine the orthogroup_id,
                            # clustering method and gene id using the format
                            # ,orthogroup_id>^^<clustering_method>^^<gene_id>.
                            combined_id = "%s^^%s^^%s" % (orthogroup_id, clustering_method, gene_id)
                            if combined_id not in adict:
                                # The value will be the dna sequence string..
                                adict[combined_id] = ""
                        else:
                            # Example line:
                            # ATGGAGAAGGACTTT
                            # Here combined_id is set because the fasta format specifies
                            # that all lines following the gene id defined in the if block
                            # above will be the sequence associated with that gene until
                            # the next gene id line is encountered.
                            sequence = adict[combined_id]
                            sequence = "%s%s" % (sequence, line)
                            adict[combined_id] = sequence
        # Populate the plant_tribes_gene and gen_scaffold_association tables
        # from the contents of aa_dict and dna_dict.
        for combined_id in sorted(dna_dict.keys()):
            # The dictionary keys combine the orthogroup_id, clustering method and
            # gene id using the format <orthogroup_id>^^<clustering_method>^^<gene_id>.
            items = combined_id.split("^^")
            orthogroup_id = items[0]
            clustering_method = items[1]
            gene_id = items[2]
            print("Populating the plant_tribes_gene and gene_scaffold_orthogroup_association tables with gene %s, scaffold %s and orthogroup %s..." % (gene_id, scaffold_id, orthogroup_id))
            # The value will be a list containing both
            # clustering_method and the dna string.
            dna_sequence = dna_dict[combined_id]
            aa_sequence = aa_dict[combined_id]
            # Get the species_code from the gene_id.
            species_code = gene_id.split("_")[1]
            # Strip the version from the species_code.
            species_code = species_code[0:5]
            # Get the species_name from self.species_ids_dict.
            species_name = self.species_ids_dict[species_code]
            # Get the plant_tribes_orthogroup primary key id  for
            # the orthogroup_id from the plant_tribes_orthogroup table.
            sql = "SELECT id FROM plant_tribes_orthogroup WHERE orthogroup_id = '%s';" % orthogroup_id
            cur = self.conn.cursor()
            cur.execute(sql)
            orthogroup_id_db = cur.fetchone()[0]
            # If the plant_tribes_gene table contains a row that has the gene_id,
            # then we'll add a row only to the gene_scaffold_orthogroup_association table.
            # Get the taxon_id  for the species_name from the plant_tribes_taxa table.
            sql = "SELECT id FROM plant_tribes_taxa WHERE species_name = '%s';" % species_name
            cur = self.conn.cursor()
            cur.execute(sql)
            taxon_id = cur.fetchone()[0]
            # If the plant_tribes_gene table contains a row that has the gene_id,
            # then we'll add a row only to the gene_scaffold_orthogroup_association table.
            sql = "SELECT id FROM plant_tribes_gene WHERE gene_id = '%s';" % gene_id
            cur = self.conn.cursor()
            cur.execute(sql)
            try:
                gene_id_db = cur.fetchone()[0]
            except:
                # Insert a row into the plant_tribes_gene table.
                args = [gene_id, taxon_id, dna_sequence, aa_sequence]
                sql = """
                    INSERT INTO plant_tribes_gene
                         VALUES (nextval('plant_tribes_gene_id_seq'), %s, %s, %s, %s)
                         RETURNING id;
                """
                cur = self._update(sql, tuple(args))
                self._flush()
                gene_id_db = cur.fetchone()[0]
            # Insert a row into the gene_scaffold_orthogroup_association table.
            # Get the scaffold_rec for the current scaffold_id and clustering_method.
            # The list is [<scaffold_id_db>, <scaffold_id>, <clustering_method>]
            for scaffold_rec in self.scaffold_recs:
                if scaffold_id in scaffold_rec and clustering_method in scaffold_rec:
                    scaffold_id_db = scaffold_rec[0]
            args = [gene_id_db, scaffold_id_db, orthogroup_id_db]
            sql = """
                INSERT INTO gene_scaffold_orthogroup_association
                     VALUES (nextval('gene_scaffold_orthogroup_association_id_seq'), %s, %s, %s);
            """
            cur = self._update(sql, tuple(args))
            self._flush()

    def process_scaffold_config_files(self):
        """
        1. Parse ~/<scaffold_id>/<scaffold_id>/.rootingOrder.config
        (e.g., ~/22Gv1.1/22Gv1.1..rootingOrder.config) to populate.
        2. Calculate the number of genes found
        for each species and add the number to self.species_genes_dict.
        3. Parse ~/<scaffold_id>/<scaffold_id>.taxaLineage.config to
        populate the plant_tribes_taxa table.
        """
        scaffold_id = os.path.basename(self.args.scaffold_id)
        print("scaffold_id: %s" % str(scaffold_id))
        file_dir = self.args.scaffold_id
        file_name = os.path.join(file_dir, '%s.rootingOrder.config' % scaffold_id)
        print("file_name: %s" % str(file_name))
        # Populate self.species_ids_dict.
        with open(file_name, "r") as fh:
            for i, line in enumerate(fh):
                line = line.strip()
                if len(line) == 0 or line.startswith("#") or line.startswith("["):
                    # Skip blank lines, comments and section headers.
                    continue
                # Example line:
                # Physcomitrella patens=Phypa
                items = line.split("=")
                self.species_ids_dict[items[1]] = items[0]
        # Get lineage information for orthogrpoup taxa.
        for scaffold_genes_dict_key in sorted(self.scaffold_genes_dict.keys()):
            # The format of key is <clustering_method>^^<orthogroup_id>.
            # For example: {"gfam^^1" : "gnl_Musac1.0_GSMUA_Achr1T11000_001"
            scaffold_genes_dict_key_items = scaffold_genes_dict_key.split("^^")
            clustering_method = scaffold_genes_dict_key_items[0]
            # Get the list of genes for the current scaffold_genes_dict_key.
            gene_list = self.scaffold_genes_dict[scaffold_genes_dict_key]
            for gene_id in gene_list:
                # Example species_code: Musac1.0, where
                # species_name is Musac and version is 1.0.
                species_code = gene_id.split("_")[1]
                # Strip the version from the species_code.
                species_code = species_code[0:5]
                # Get the species_name from self.species_ids_dict.
                species_name = self.species_ids_dict[species_code]
                # Create a key for self.species_genes_dict, with the format:
                # <clustering_method>^^<species_code>
                species_genes_dict_key = "%s^^%s" % (clustering_method, species_code)
                # Add an entry to self.species_genes_dict, where the value
                # is a list containing species_name and num_genes.
                if species_genes_dict_key in self.species_genes_dict:
                    tup = self.species_genes_dict[species_genes_dict_key]
                    tup[1] += 1
                    self.species_genes_dict[species_genes_dict_key] = tup
                else:
                    self.species_genes_dict[species_genes_dict_key] = [species_name, 1]
        # Populate the plant_tribes_taxa table.
        file_name = os.path.join(file_dir, '%s.taxaLineage.config' % scaffold_id)
        print("file_name: %s" % str(file_name))
        with open(file_name, "r") as fh:
            for i, line in enumerate(fh):
                line = line.strip()
                if len(line) == 0 or line.startswith("#") or line.startswith("Species"):
                    # Skip blank lines, comments and section headers.
                    continue
                # Example line: Populus trichocarpa\tSalicaceae\tMalpighiales\tRosids\tCore Eudicots
                items = line.split("\t")
                species_name = items[0]
                print("Calculating the number of genes for species_name: %s" % str(species_name))
                for species_genes_dict_key in sorted(self.species_genes_dict.keys()):
                    # The format of species_genes_dict_key is <clustering_method>^^<species_code>.
                    species_genes_dict_key_items = species_genes_dict_key.split("^^")
                    clustering_method = species_genes_dict_key_items[0]
                    species_code = species_genes_dict_key_items[1]
                    # Get the scaffold_rec for the current scaffold_id and clustering_method.
                    # The list is [<scaffold_id_db>, <scaffold_id>, <clustering_method>]
                    for scaffold_rec in self.scaffold_recs:
                        if scaffold_id in scaffold_rec and clustering_method in scaffold_rec:
                            scaffold_id_db = scaffold_rec[0]
                    # The value is a list containing species_name and num_genes.
                    val = self.species_genes_dict[species_genes_dict_key]
                    if species_name == val[0]:
                        num_genes = val[1]
                    else:
                        num_genes = 0
                    # Insert a row in to the plant_tribes_scaffold table.
                    args = [species_name, scaffold_id_db, num_genes, items[1], items[2], items[3], items[4]]
                    sql = """
                        INSERT INTO plant_tribes_taxa
                             VALUES (nextval('plant_tribes_taxa_id_seq'), %s, %s, %s, %s, %s, %s, %s);
                    """
                    self._update(sql, tuple(args))
                    self._flush()

    def process_annot_dir(self):
        """
        First, parse all of the *.min_evalue.summary files in the
        ~/<scaffold_id>/annot directory (e.g., ~/22Gv1.1/annot) to populate
        both the plant_tribes_scaffold and the plant_tribes_orthogroup tables.
        Next, parse all of the *.list files in the same directory to populate
        self.scaffold_genes_dict.
        """
        scaffold_id = os.path.basename(self.args.scaffold_id)
        file_dir = os.path.join(self.args.scaffold_id, 'annot')
        # The scaffol naming convention must follow this pattern:
        # <integer1>Gv<integer2>.<integer3>
        # where integer 1 is the number of genomes in the scaffold_id.  For example:
        # 22Gv1.1 -> 22 genomes
        # 12Gv1.0 -> 12 genomes
        # 26Gv2.0 -> 26 genomes, etc.
        num_genomes = int(scaffold_id.split("Gv")[0])
        super_ortho_start_index = num_genomes + 1
        for file_name in glob.glob(os.path.join(file_dir, "*min_evalue.summary")):
            items = os.path.basename(file_name).split(".")
            clustering_method = items[0]
            # Save all clustering methods for later processing.
            if clustering_method not in self.clustering_methods:
                self.clustering_methods.append(clustering_method)
            # Insert a row in to the plant_tribes_scaffold table.
            print("Inserting a row into the plant_tribes_scaffold table for scaffold %s and clustering method %s..." % (scaffold_id, clustering_method))
            args = [scaffold_id, clustering_method]
            sql = """
                INSERT INTO plant_tribes_scaffold
                     VALUES (nextval('plant_tribes_scaffold_id_seq'), %s, %s)
                     RETURNING id;
            """
            cur = self._update(sql, tuple(args))
            self._flush()
            scaffold_id_db = cur.fetchone()[0]
            self.scaffold_recs.append([scaffold_id_db, scaffold_id, clustering_method])
            with open(file_name, "r") as fh:
                for i, line in enumerate(fh):
                    if i == 0:
                        # Skip first line.
                        continue
                    num_genes = 0
                    num_species = 0
                    items = line.split("\t")
                    orthogroup_id = int(items[0])
                    # Zero based items 1 to num_genomes consists of the
                    # number of species classified in the orthogroup (i.e.,
                    # species with at least 1 gene in the orthogroup).
                    for j in range(1, num_genomes):
                        j_int = int(items[j])
                        if j_int > 0:
                            # The  species has at least 1 gene
                            num_species += 1
                            num_genes += j_int
                    # Insert a row into the plant_tribes_orthogroup table.
                    print("Inserting a row into the plant_tribes_orthogroup table...")
                    args = [orthogroup_id, scaffold_id_db, num_species, num_genes]
                    for k in range(super_ortho_start_index, len(items)):
                        args.append('%s' % str(items[k]))
                    sql = """
                        INSERT INTO plant_tribes_orthogroup
                             VALUES (nextval('plant_tribes_orthogroup_id_seq'), %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);
                    """
                    cur = self._update(sql, tuple(args))
                    self._flush()
        for file_name in glob.glob(os.path.join(file_dir, "*list")):
            items = os.path.basename(file_name).split(".")
            clustering_method = items[0]
            with open(file_name, "r") as fh:
                for i, line in enumerate(fh):
                    items = line.split("\t")
                    # The key will be a combination of clustering_method and
                    # orthogroup_id separated by "^^" for easy splitting later.
                    key = "%s^^%s" % (clustering_method, items[0])
                    # The value is the gen_id with all white space replaced by "_".
                    val = items[1].replace("|", "_")
                    if key in self.scaffold_genes_dict:
                        self.scaffold_genes_dict[key].append(val)
                    else:
                        self.scaffold_genes_dict[key] = [val]

    def _flush(self):
        if self.args.dry_run:
            self.conn.rollback()
            print("--dry-run specified, all changes rolled back")
        else:
            self.conn.commit()

    def _shutdown(self):
        self.conn.close()


if __name__ == '__main__':
    add_scaffold = AddScaffold()
    add_scaffold._run()
    add_scaffold._shutdown()

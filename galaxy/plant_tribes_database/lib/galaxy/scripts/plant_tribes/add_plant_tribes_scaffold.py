#!/usr/bin/env python
"""
add_plant_tribes_scaffold.py - A script for adding a scaffold to the Galaxy PlantTribes
database efficiently by bypassing the Galaxy model and operating directly on the database.
PostgreSQL 9.1 or greater is required.
"""
from __future__ import print_function

import argparse
import datetime
import glob
import inspect
import logging
import os
import shutil
import sys

import psycopg2
from sqlalchemy.engine.url import make_url

galaxy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
sys.path.insert(1, os.path.join(galaxy_root, 'lib'))

import galaxy.config
from galaxy.util.script import app_properties_from_args, populate_config_args

log = logging.getLogger()


class AddScaffold(object):
    def __init__(self):
        self.args = None
        self.config = None
        self.conn = None
        self.logs = {}
        self.__parse_args()
        self.__setup_logging()
        self.__load_config()
        self.__connect_db()

    def __parse_args(self):
        parser = argparse.ArgumentParser()
        populate_config_args(parser)
        parser.add_argument('--debug', action='store_true', dest='debug', help='Enable debug logging', default=False)
        parser.add_argument('--dry_run', action='store_true', dest='dry_run', help="Dry run (rollback all transactions)", default=False)
        parser.add_argument('--log_dir', dest='log_dir', help='Log file directory', default=os.path.join(galaxy_root, 'scripts', 'plant_tribes'))
        parser.add_argument('--scaffold', dest='scaffold', help='Full path to PlantTribes scaffold directory', default=os.path.join(galaxy_root, 'plant_tribes', 'scaffolds', '22Gv1.1'))
        parser.add_argument('--taxa_lineage', dest='taxa_lineage', help='Scaffold taxa lineage', default=os.path.join(galaxy_root, 'plant_tribes', 'scaffolds', '22Gv1.1', '22Gv1.1taxaLineage'))
        self.args = parser.parse_args()

    def __setup_logging(self):
        format = "%(funcName)s %(levelname)s %(asctime)s %(message)s"
        if self.args.debug:
            logging.basicConfig(level=logging.DEBUG, format=format)
        else:
            logging.basicConfig(level=logging.INFO, format=format)

    def __load_config(self):
        app_properties = app_properties_from_args(self.args)
        self.config = galaxy.config.Configuration(**app_properties)

    def __connect_db(self):
        url = make_url(galaxy.config.get_plant_tribes_database_url(self.config))
        log.info('Connecting to database with URL: %s' % url)
        args = url.translate_connect_args(username='user')
        args.update(url.query)

        assert url.get_dialect().name == 'postgresql', 'This script can only be used with PostgreSQL.'

        self.conn = psycopg2.connect(**args)

    def _open_logfile(self):
        action_name = inspect.stack()[1][3]
        logname = os.path.join(self.args.log_dir, action_name + '.log')

        if self.args.dry_run:
            log.debug('--dry-run specified, logging changes to stdout instead of log file: %s' % logname)
            self.logs[action_name] = sys.stdout
        else:
            log.debug('Opening log file: %s' % logname)
            self.logs[action_name] = open(logname, 'a')

        message = '==== Log opened: %s ' % datetime.datetime.now().isoformat()
        self.logs[action_name].write(message.ljust(72, '='))
        self.logs[action_name].write('\n')

    def _log(self, message, action_name=None):
        if action_name is None:
            action_name = inspect.stack()[1][3]
        if not message.endswith('\n'):
            message += '\n'
        self.logs[action_name].write(message)

    def _close_logfile(self):
        action_name = inspect.stack()[1][3]

        message = '==== Log closed: %s ' % datetime.datetime.now().isoformat()
        self.logs[action_name].write(message.ljust(72, '='))
        self.logs[action_name].write('\n')

        if self.args.dry_run:
            log.debug('--dry-run specified, changes were logged to stdout insted of log file')
        else:
            log.debug('Closing log file: %s' % self.logs[action_name].name)
            self.logs[action_name].close()

        del self.logs[action_name]

    def _run(self):
        self.process_scaffold_annot_dir()

    def _update(self, sql, args):
        if args is not None:
            log.debug('SQL is: %s' % sql % args)
            print('SQL is: %s' % sql % args)
        else:
            log.debug('SQL is: %s' % sql)
        log.debug('SQL is: %s' % sql)
        log.debug('args is: %s' % str(args))
        cur = self.conn.cursor()
        log.info('Executing SQL')
        cur.execute(sql, args)
        log.info('Database status: %s' % cur.statusmessage)

        return cur

    def process_scaffold_annot_dir(self):
        """
        Process the files in scaffold annot directory ~/22Gv1.1/annot.
        """
        annot_dir = os.path.join(self.args.scaffold, 'annot')
        print("annot_dir: %s" % str(annot_dir))
        scaffold = os.path.basename(self.args.scaffold)
        print("scaffold: %s" % str(scaffold))
        # The scaffol naming convention must follow this pattern:
        # <integer1>Gv<integer2>.<integer3>
        # where integer 1 is the number of genomes in the scaffold.  For example:
        # 22Gv1.1 -> 22 genomes
        # 12Gv1.0 -> 12 genomes
        # 26Gv2.0 -> 26 genomes, etc.
        num_genomes = int(scaffold.split("Gv")[0])
        print("num_genomes: %s" % str(num_genomes))
        super_ortho_start_index = num_genomes + 1
        print("super_ortho_start_index: %s" % str(super_ortho_start_index))
        for file_name in glob.glob(os.path.join(annot_dir, "*min_evalue.summary")):
            print("filename: %s" % str(file_name))
            items = os.path.basename(file_name).split(".")
            print("items: %s" % str(items))
            clustering_method = items[0]
            with open(file_name, "r") as fh:
                for i, line in enumerate(fh):
                    if i== 0:
                        # Skip first line.
                        continue
                    num_genes = 0
                    num_species = 0
                    items = line.split("\t")
                    ortho_id = int(items[0])
                    # Zero based items 1 to num_genomes consists of the
                    # number of species classified in the orthogroup (i.e.,
                    # species with at least 1 gene in the orthogroup).
                    for j in range(1, num_genomes):
                        j_int = int(items[j])
                        if j_int > 0:
                            # The  species has at least 1 gene
                            num_species += 1
                            num_genes += j_int
                    # TODO: fix this to insert into orthogroups table.
                    """
                    Column("ortho_id", Integer, primary_key=True, nullable=False),
                    Column("scaffold_id", TrimmedString(10), ForeignKey("plant_tribes_scaffold.scaffold_id"), index=True, nullable=False),
                    Column("clustering_method", TrimmedString(30), nullable=False),
                    Column("num_taxa", Integer, nullable=False),
                    Column("num_genes", Integer, nullable=False),
                    Column("super_ortho_1_2", TrimmedString(10), nullable=False),
                    Column("super_ortho_1_5", TrimmedString(10), nullable=False),
                    Column("super_ortho_1_8", TrimmedString(10), nullable=False),
                    Column("super_ortho_2_0", TrimmedString(10), nullable=False),
                    Column("super_ortho_2_5", TrimmedString(10), nullable=False),
                    Column("super_ortho_3_0", TrimmedString(10), nullable=False),
                    Column("super_ortho_3_5", TrimmedString(10), nullable=False),
                    Column("super_ortho_4_0", TrimmedString(10), nullable=False),
                    Column("super_ortho_4_5", TrimmedString(10), nullable=False),
                    Column("super_ortho_5_0", TrimmedString(10), nullable=False),
                    Column("ahdr_description", TEXT, index=True, nullable=False),
                    Column("tair_description", TEXT, index=True, nullable=False),
                    Column("pfam_description", TEXT, index=True, nullable=False),
                    Column("interproscan__description", TEXT, index=True, nullable=False),
                    Column("gene_ontology_molecular_funcion_description", TEXT, index=True, nullable=False),
                    Column("gene_ontology_biological_process_description", TEXT, index=True, nullable=False),
                    Column("gene_ontology_celular_component_description", TEXT, index=True, nullable=False))

                    Eric's Perl code:
                    # output for orthogroups database table file
                    print OUT "$fields[0]\t$scaffold\t$clustering_method\t$num_of_species\t$num_of_genes";
                    my $super_ortho_start = $genomes + 1;
                    for ($super_ortho_start..$#fields){
                        print OUT "\t$fields[$_]";
                    }
                    print OUT "\n";
                    """
                    args = [ortho_id, scaffold, clustering_method, num_species, num_genes]
                    #log.debug('args 1: %s' % str(args))
                    print('args 1: %s' % str(args))
                    #print("items: %s" % str(items))
                    print("items[super_ortho_start_index]: %s" % str(items[super_ortho_start_index]))
                    print("len(items): %s" % str(len(items)))
                    print("len(items)-1: %s" % str(len(items)-1))
                    print("items[len(items)-1]: %s" % str(items[len(items)-1]))
                    for k in range(super_ortho_start_index, len(items)-1):
                        print("k: %s" % str(k))
                        print("items[k]: %s" % str(items[k]))
                        log.debug('k: %s' % str(k))
                        log.debug('items[k]: %s' % str(items[k]))
                        args.append('%s' % str(items[k]))
                    log.debug('args 2: %s' % str(args))
                    print('args 2: %s' % str(args))
                    print('len(args) 2: %s' % str(len(args)))
                    sql = """
                        INSERT INTO plant_tribes_orthogroup
                             VALUES (%d, '%s', '%s', %d, %d, '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s');
                    """
                    cur = self._update(sql, tuple(args))
                    self._flush()

    def _flush(self):
        if self.args.dry_run:
            self.conn.rollback()
            log.info("--dry-run specified, all changes rolled back")
        else:
            self.conn.commit()
            log.info("All changes committed")

    def _shutdown(self):
        self.conn.close()
        for handle in self.logs.values():
            message = '==== Log closed at shutdown: %s ' % datetime.datetime.now().isoformat()
            handle.write(message.ljust(72, '='))
            handle.write('\n')
            handle.close()


if __name__ == '__main__':
    add_scaffold = AddScaffold()
    try:
        add_scaffold._run()
    except Exception:
        log.exception('Caught exception in run sequence:')
    add_scaffold._shutdown()

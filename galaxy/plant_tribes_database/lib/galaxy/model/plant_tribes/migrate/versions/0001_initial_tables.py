import datetime
import logging

from sqlalchemy import Column, DateTime, ForeignKey, Integer, MetaData, Table, TEXT

# Need our custom types, but don't import anything else from model
from galaxy.model.custom_types import TrimmedString

now = datetime.datetime.utcnow
log = logging.getLogger(__name__)
metadata = MetaData()

# Tables as of changeset 1464:c7acaa1bb88f
PlantTribesScaffold_table = Table("plant_tribes_scaffold", metadata,
                                  Column("id", Integer, primary_key=True),
                                  Column("create_time", DateTime, default=now),
                                  Column("update_time", DateTime, default=now, onupdate=now),
                                  Column("description", TEXT, nullable=False),
                                  Column("genes", Integer, nullable=False),
                                  Column("orthogroups", Integer))

PlantTribesTaxa_table = Table("plant_tribes_taxa", metadata,
                              Column("id", Integer, primary_key=True),
                              Column("create_time", DateTime, default=now),
                              Column("update_time", DateTime, default=now, onupdate=now),
                              Column("genes", Integer, nullable=False),
                              Column("species_family", TrimmedString(50), nullable=False),
                              Column("species_order", TrimmedString(50), nullable=False),
                              Column("species_group", TrimmedString(50), nullable=False),
                              Column("species_clade", TrimmedString(50), nullable=False))

PlantTribesOrthogroup_table = Table("plant_tribes_orthogroup", metadata,
                                    Column("id", Integer, primary_key=True),
                                    Column("create_time", DateTime, default=now),
                                    Column("update_time", DateTime, default=now, onupdate=now),
                                    Column("taxa", Integer, nullable=False),
                                    Column("genes", Integer, nullable=False),
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
                                    Column("ahdr_desc", TEXT, nullable=False),
                                    Column("tair_desc", TEXT, nullable=False),
                                    Column("pfam_desc", TEXT, nullable=False),
                                    Column("iprscan_desc", TEXT, nullable=False),
                                    Column("go_mf_desc", TEXT, nullable=False),
                                    Column("go_bp_desc", TEXT, nullable=False),
                                    Column("go_cc_desc", TEXT, nullable=False))

PlantTribesGene_table = Table("plant_tribes_gene", metadata,
                              Column("id", Integer, primary_key=True),
                              Column("create_time", DateTime, default=now),
                              Column("update_time", DateTime, default=now, onupdate=now),
                              Column("dna_sequence", TEXT, nullable=False),
                              Column("aa_sequence", TEXT, nullable=False))

ScaffoldAssociation_table = Table("scaffold_association", metadata,
                                  Column("id", Integer, primary_key=True),
                                  Column("scaffold", TrimmedString(10), index=True, nullable=False),
                                  Column("scaffold_id", Integer, ForeignKey("plant_tribes_scaffold.id"), index=True, nullable=False),
                                  Column("taxa_id", Integer, ForeignKey("plant_tribes_taxa.id"), index=True, nullable=False),
                                  Column("orthogroup_id", Integer, ForeignKey("plant_tribes_orthogroup.id"), index=True, nullable=False),
                                  Column("gene_id", Integer, ForeignKey("plant_tribes_gene.id"), index=True, nullable=False))

SpeciesAssociation_table = Table("species_association", metadata,
                                 Column("id", Integer, primary_key=True),
                                 Column("species", TrimmedString(50), index=True, nullable=False),
                                 Column("taxa_id", Integer, ForeignKey("plant_tribes_taxa.id"), index=True, nullable=False),
                                 Column("gene_id", Integer, ForeignKey("plant_tribes_gene.id"), index=True, nullable=False))


def upgrade(migrate_engine):
    metadata.bind = migrate_engine
    metadata.create_all()

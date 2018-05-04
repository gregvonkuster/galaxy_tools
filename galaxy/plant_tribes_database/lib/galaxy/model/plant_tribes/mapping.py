from sqlalchemy import (
    Column,
    DateTime,
    ForeignKey,
    Integer,
    MetaData,
    Table,
    TEXT
)
from sqlalchemy.orm import (
    mapper,
    relation
)

from galaxy.model import plant_tribes as plant_tribes_model
from galaxy.model.base import ModelMapping
from galaxy.model.custom_types import TrimmedString
from galaxy.model.orm.engine_factory import build_engine
from galaxy.model.orm.now import now

metadata = MetaData()

plant_tribes_model.PlantTribesScaffold.table = Table("plant_tribes_scaffold", metadata,
                                                     Column("id", Integer, primary_key=True),
                                                     Column("create_time", DateTime, default=now),
                                                     Column("update_time", DateTime, default=now, onupdate=now),
                                                     Column("description", TEXT, nullable=False),
                                                     Column("genes", Integer, nullable=False),
                                                     Column("orthogroups", Integer))

plant_tribes_model.PlantTribesTaxa.table = Table('plant_tribes_taxa', metadata,
                                                 Column("id", Integer, primary_key=True),
                                                 Column("create_time", DateTime, default=now),
                                                 Column("update_time", DateTime, default=now, onupdate=now),
                                                 Column("genes", Integer, nullable=False),
                                                 Column("species_family", TrimmedString(50), nullable=False),
                                                 Column("species_order", TrimmedString(50), nullable=False),
                                                 Column("species_group", TrimmedString(50), nullable=False),
                                                 Column("species_clade", TrimmedString(50), nullable=False))

plant_tribes_model.PlantTribesOrthogroup.table = Table("plant_tribes_orthogroup", metadata,
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

plant_tribes_model.PlantTribesGene.table = Table("plant_tribes_gene", metadata,
                                                 Column("id", Integer, primary_key=True),
                                                 Column("create_time", DateTime, default=now),
                                                 Column("update_time", DateTime, default=now, onupdate=now),
                                                 Column("dna_sequence", TEXT, nullable=False),
                                                 Column("aa_sequence", TEXT, nullable=False))

plant_tribes_model.ScaffoldAssociation.table = Table("scaffold_association", metadata,
                                                     Column("id", Integer, primary_key=True),
                                                     Column("scaffold", TrimmedString(10), index=True, nullable=False),
                                                     Column("scaffold_id", Integer, ForeignKey("plant_tribes_scaffold.id"), index=True, nullable=False),
                                                     Column("taxa_id", Integer, ForeignKey("plant_tribes_taxa.id"), index=True, nullable=False),
                                                     Column("orthogroup_id", Integer, ForeignKey("plant_tribes_orthogroup.id"), index=True, nullable=False),
                                                     Column("gene_id", Integer, ForeignKey("plant_tribes_gene.id"), index=True, nullable=False))

plant_tribes_model.SpeciesAssociation.table = Table("species_association", metadata,
                                                    Column("id", Integer, primary_key=True),
                                                    Column("species", TrimmedString(50), index=True, nullable=False),
                                                    Column("taxa_id", Integer, ForeignKey("plant_tribes_taxa.id"), index=True, nullable=False),
                                                    Column("gene_id", Integer, ForeignKey("plant_tribes_gene.id"), index=True, nullable=False))


mapper(plant_tribes_model.PlantTribesScaffold, plant_tribes_model.PlantTribesScaffold.table,
       properties=dict(scaffolds=relation(plant_tribes_model.ScaffoldAssociation,
                                          primaryjoin=(plant_tribes_model.PlantTribesScaffold.table.c.id == plant_tribes_model.ScaffoldAssociation.table.c.scaffold_id))))

mapper(plant_tribes_model.PlantTribesTaxa, plant_tribes_model.PlantTribesTaxa.table,
       properties=dict(scaffolds=relation(plant_tribes_model.ScaffoldAssociation,
                                          primaryjoin=(plant_tribes_model.PlantTribesTaxa.table.c.id == plant_tribes_model.ScaffoldAssociation.table.c.taxa_id)),
                       species=relation(plant_tribes_model.SpeciesAssociation,
                                        primaryjoin=(plant_tribes_model.PlantTribesTaxa.table.c.id == plant_tribes_model.SpeciesAssociation.table.c.taxa_id))))

mapper(plant_tribes_model.PlantTribesOrthogroup, plant_tribes_model.PlantTribesOrthogroup.table,
       properties=dict(scaffolds=relation(plant_tribes_model.ScaffoldAssociation,
                                          primaryjoin=(plant_tribes_model.PlantTribesOrthogroup.table.c.id == plant_tribes_model.ScaffoldAssociation.table.c.orthogroup_id))))

mapper(plant_tribes_model.PlantTribesGene, plant_tribes_model.PlantTribesGene.table,
       properties=dict(scaffolds=relation(plant_tribes_model.ScaffoldAssociation,
                                          primaryjoin=(plant_tribes_model.PlantTribesGene.table.c.id == plant_tribes_model.ScaffoldAssociation.table.c.gene_id)),
                       species=relation(plant_tribes_model.SpeciesAssociation,
                                        primaryjoin=(plant_tribes_model.PlantTribesGene.table.c.id == plant_tribes_model.SpeciesAssociation.table.c.gene_id))))

mapper(plant_tribes_model.ScaffoldAssociation, plant_tribes_model.ScaffoldAssociation.table,
       properties=dict(plant_tribes_scaffolds=relation(plant_tribes_model.PlantTribesScaffold,
                                                       primaryjoin=(plant_tribes_model.ScaffoldAssociation.table.c.scaffold_id == plant_tribes_model.PlantTribesScaffold.table.c.id)),
                       plant_tribes_taxa=relation(plant_tribes_model.PlantTribesTaxa,
                                                  primaryjoin=(plant_tribes_model.ScaffoldAssociation.table.c.taxa_id == plant_tribes_model.PlantTribesTaxa.table.c.id)),
                       plant_tribes_orthogroups=relation(plant_tribes_model.PlantTribesOrthogroup,
                                                         primaryjoin=(plant_tribes_model.ScaffoldAssociation.table.c.orthogroup_id == plant_tribes_model.PlantTribesOrthogroup.table.c.id)),
                       plant_tribes_genes=relation(plant_tribes_model.PlantTribesGene,
                                                   primaryjoin=(plant_tribes_model.ScaffoldAssociation.table.c.gene_id == plant_tribes_model.PlantTribesGene.table.c.id))))

mapper(plant_tribes_model.SpeciesAssociation, plant_tribes_model.SpeciesAssociation.table,
       properties=dict(plant_tribes_taxa=relation(plant_tribes_model.PlantTribesTaxa,
                                                  primaryjoin=(plant_tribes_model.SpeciesAssociation.table.c.taxa_id == plant_tribes_model.PlantTribesTaxa.table.c.id)),
                       plant_tribes_genea=relation(plant_tribes_model.PlantTribesGene,
                                                   primaryjoin=(plant_tribes_model.SpeciesAssociation.table.c.gene_id == plant_tribes_model.PlantTribesGene.table.c.id))))


def init(url, engine_options={}, create_tables=False):
    """Connect mappings to the database"""
    # Load the appropriate db module
    engine = build_engine(url, engine_options)
    # Connect the metadata to the database.
    metadata.bind = engine
    result = ModelMapping([plant_tribes_model], engine=engine)
    # Create tables if needed
    if create_tables:
        metadata.create_all()
        # metadata.engine.commit()
    result.create_tables = create_tables
    # load local galaxy security policy
    return result

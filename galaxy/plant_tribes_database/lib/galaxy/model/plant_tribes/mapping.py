from sqlalchemy import (
    Column,
    ForeignKey,
    Integer,
    JSON,
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

metadata = MetaData()

plant_tribes_model.PlantTribesScaffold.table = Table("plant_tribes_scaffold", metadata,
                                                     Column("id", Integer, primary_key=True),
                                                     Column("scaffold_id", TrimmedString(10), index=True, nullable=False),
                                                     Column("clustering_method", TrimmedString(30), index=True, nullable=False))

plant_tribes_model.PlantTribesTaxon.table = Table('plant_tribes_taxon', metadata,
                                                  Column("id", Integer, primary_key=True),
                                                  Column("species_name", TrimmedString(50), index=True, nullable=False),
                                                  Column("scaffold_id", Integer, ForeignKey("plant_tribes_scaffold.id"), index=True, nullable=False),
                                                  Column("num_genes", Integer, nullable=False),
                                                  Column("species_family", TrimmedString(50), nullable=False),
                                                  Column("species_order", TrimmedString(50), nullable=False),
                                                  Column("species_group", TrimmedString(50), nullable=False),
                                                  Column("species_clade", TrimmedString(50), nullable=False))

plant_tribes_model.PlantTribesOrthogroup.table = Table("plant_tribes_orthogroup", metadata,
                                                       Column("id", Integer, primary_key=True),
                                                       Column("orthogroup_id", Integer, index=True, nullable=False),
                                                       Column("scaffold_id", Integer, ForeignKey("plant_tribes_scaffold.id"), index=True, nullable=False),
                                                       Column("num_species", Integer, nullable=False),
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
                                                       Column("ahdr_description", JSON, nullable=False),
                                                       Column("tair_description", JSON, nullable=False),
                                                       Column("pfam_description", JSON, nullable=False),
                                                       Column("interproscan_description", JSON, nullable=False),
                                                       Column("molecular_function", JSON, nullable=False),
                                                       Column("biological_process", JSON, nullable=False),
                                                       Column("cellular_component", JSON, nullable=False))

plant_tribes_model.PlantTribesGene.table = Table("plant_tribes_gene", metadata,
                                                 Column("id", Integer, primary_key=True),
                                                 Column("gene_id", TrimmedString(100), index=True, nullable=False),
                                                 Column("dna_sequence", TEXT, nullable=False),
                                                 Column("aa_sequence", TEXT, nullable=False))

plant_tribes_model.GeneScaffoldOrthogroupTaxonAssociation.table = Table("gene_scaffold_orthogroup_taxon_association", metadata,
                                                                        Column("id", Integer, primary_key=True),
                                                                        Column("gene_id",Integer, ForeignKey("plant_tribes_gene.id"), index=True, nullable=False),
                                                                        Column("scaffold_id", Integer, ForeignKey("plant_tribes_scaffold.id"), index=True, nullable=False),
                                                                        Column("orthogroup_id", Integer, ForeignKey("plant_tribes_orthogroup.id"), index=True, nullable=False),
                                                                        Column("taxon_id", Integer, ForeignKey("plant_tribes_taxon.id"), index=True, nullable=False))

mapper(plant_tribes_model.PlantTribesScaffold, plant_tribes_model.PlantTribesScaffold.table,
       properties=dict(orthogroup=relation(plant_tribes_model.PlantTribesOrthogroup,
                                           backref="scaffold",
                                           primaryjoin=(plant_tribes_model.PlantTribesScaffold.table.c.id == plant_tribes_model.PlantTribesOrthogroup.table.c.scaffold_id))))

mapper(plant_tribes_model.PlantTribesTaxon, plant_tribes_model.PlantTribesTaxon.table,
       properties=dict(scaffold=relation(plant_tribes_model.PlantTribesScaffold,
                                         backref="taxa",
                                         primaryjoin=(plant_tribes_model.PlantTribesTaxon.table.c.scaffold_id == plant_tribes_model.PlantTribesScaffold.table.c.id))))

mapper(plant_tribes_model.PlantTribesOrthogroup, plant_tribes_model.PlantTribesOrthogroup.table)

mapper(plant_tribes_model.PlantTribesGene, plant_tribes_model.PlantTribesGene.table)

mapper(plant_tribes_model.GeneScaffoldOrthogroupTaxonAssociation, plant_tribes_model.GeneScaffoldOrthogroupTaxonAssociation.table,
       properties=dict(gene=relation(plant_tribes_model.PlantTribesGene),
                       scaffold=relation(plant_tribes_model.PlantTribesScaffold,
                                         lazy=False,
                                         backref="gene"),
                       orthogroup=relation(plant_tribes_model.PlantTribesOrthogroup,
                                           lazy=False,
                                           backref="gene"),
                       taxon=relation(plant_tribes_model.PlantTribesTaxon,
                                      lazy=False,
                                      backref="gene")))


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

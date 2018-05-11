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
                                                     Column("scaffold_id", TrimmedString(10), primary_key=True),
                                                     Column("description", TEXT, nullable=False),
                                                     Column("num_genes", Integer, nullable=False),
                                                     Column("num_orthogroups", Integer))

plant_tribes_model.PlantTribesTaxa.table = Table('plant_tribes_taxa', metadata,
                                                 Column("scaffold_id", TrimmedString(10), ForeignKey("plant_tribes_scaffold.scaffold_id"), primary_key=True, nullable=False),
                                                 Column("num_genes", Integer, nullable=False),
                                                 Column("species", TrimmedString(50), index=True, nullable=False),
                                                 Column("species_family", TrimmedString(50), nullable=False),
                                                 Column("species_order", TrimmedString(50), nullable=False),
                                                 Column("species_group", TrimmedString(50), nullable=False),
                                                 Column("species_clade", TrimmedString(50), nullable=False))

plant_tribes_model.PlantTribesOrthogroup.table = Table("plant_tribes_orthogroup", metadata,
                                                       Column("ortho_id", Integer, primary_key=True, nullable=False),
                                                       Column("scaffold_id", TrimmedString(10), ForeignKey("plant_tribes_scaffold.scaffold_id"), index=True, nullable=False),
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

plant_tribes_model.PlantTribesGene.table = Table("plant_tribes_gene", metadata,
                                                 Column("gene_id", Integer, primary_key=True, nullable=False),
                                                 Column("species", TrimmedString(50), nullable=False),
                                                 Column("ortho_id", Integer, ForeignKey("plant_tribes_orthogroup.ortho_id"), index=True, nullable=False),
                                                 Column("scaffold_id", TrimmedString(10), ForeignKey("plant_tribes_scaffold.scaffold_id"), index=True, nullable=False),
                                                 Column("dna_sequence", TEXT, nullable=False),
                                                 Column("aa_sequence", TEXT, nullable=False))

mapper(plant_tribes_model.PlantTribesScaffold, plant_tribes_model.PlantTribesScaffold.table,
       properties=dict(taxa=relation(plant_tribes_model.PlantTribesTaxa,
                                     primaryjoin=(plant_tribes_model.PlantTribesScaffold.table.c.scaffold_id == plant_tribes_model.PlantTribesTaxa.table.c.scaffold_id)),
                       orthogroup=relation(plant_tribes_model.PlantTribesOrthogroup,
                                     primaryjoin=(plant_tribes_model.PlantTribesScaffold.table.c.scaffold_id == plant_tribes_model.PlantTribesOrthogroup.table.c.scaffold_id))))

mapper(plant_tribes_model.PlantTribesTaxa, plant_tribes_model.PlantTribesTaxa.table,
       properties=dict(scaffold=relation(plant_tribes_model.PlantTribesScaffold,
                                         primaryjoin=(plant_tribes_model.PlantTribesTaxa.table.c.scaffold_id == plant_tribes_model.PlantTribesScaffold.table.c.scaffold_id))))

mapper(plant_tribes_model.PlantTribesOrthogroup, plant_tribes_model.PlantTribesOrthogroup.table,
       properties=dict(scaffold=relation(plant_tribes_model.PlantTribesScaffold,
                                         primaryjoin=(plant_tribes_model.PlantTribesOrthogroup.table.c.scaffold_id == plant_tribes_model.PlantTribesScaffold.table.c.scaffold_id))))

mapper(plant_tribes_model.PlantTribesGene, plant_tribes_model.PlantTribesGene.table,
       properties=dict(orthogroup=relation(plant_tribes_model.PlantTribesOrthogroup,
                                           primaryjoin=(plant_tribes_model.PlantTribesGene.table.c.ortho_id == plant_tribes_model.PlantTribesOrthogroup.table.c.ortho_id)),
                       scaffold=relation(plant_tribes_model.PlantTribesScaffold,
                                         primaryjoin=(plant_tribes_model.PlantTribesGene.table.c.scaffold_id == plant_tribes_model.PlantTribesScaffold.table.c.scaffold_id))))


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

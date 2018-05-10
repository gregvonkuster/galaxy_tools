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
                                                     Column("scaffold", TrimmedString(10), index=True, nullable=False),
                                                     Column("description", TEXT, nullable=False),
                                                     Column("num_genes", Integer, nullable=False),
                                                     Column("num_orthogroups", Integer))

plant_tribes_model.PlantTribesTaxa.table = Table('plant_tribes_taxa', metadata,
                                                 Column("id", Integer, primary_key=True),
                                                 Column("create_time", DateTime, default=now),
                                                 Column("update_time", DateTime, default=now, onupdate=now),
                                                 Column("num_genes", Integer, nullable=False),
                                                 Column("species", TrimmedString(50), index=True, nullable=False),
                                                 Column("species_family", TrimmedString(50), nullable=False),
                                                 Column("species_order", TrimmedString(50), nullable=False),
                                                 Column("species_group", TrimmedString(50), nullable=False),
                                                 Column("species_clade", TrimmedString(50), nullable=False))

plant_tribes_model.PlantTribesOrthogroup.table = Table("plant_tribes_orthogroup", metadata,
                                                       Column("id", Integer, primary_key=True),
                                                       Column("create_time", DateTime, default=now),
                                                       Column("update_time", DateTime, default=now, onupdate=now),
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
                                                 Column("id", Integer, primary_key=True),
                                                 Column("create_time", DateTime, default=now),
                                                 Column("update_time", DateTime, default=now, onupdate=now),
                                                 Column("dna_sequence", TEXT, nullable=False),
                                                 Column("aa_sequence", TEXT, nullable=False))

plant_tribes_model.TaxaScaffoldAssociation.table = Table("taxa_scaffold_association", metadata,
                                                     Column("id", Integer, primary_key=True),
                                                     Column("taxa_id", Integer, ForeignKey("plant_tribes_taxa.id"), index=True, nullable=False),
                                                     Column("scaffold_id", Integer, ForeignKey("plant_tribes_scaffold.id"), index=True, nullable=False))

plant_tribes_model.TaxaGeneAssociation.table = Table("taxa_gene_association", metadata,
                                                    Column("id", Integer, primary_key=True),
                                                    Column("taxa_id", Integer, ForeignKey("plant_tribes_taxa.id"), index=True, nullable=False),
                                                    Column("gene_id", Integer, ForeignKey("plant_tribes_gene.id"), index=True, nullable=False))

plant_tribes_model.ScaffoldAssociation.table = Table("scaffold_taxa_orthogroup_gene_association", metadata,
                                                     Column("id", Integer, primary_key=True),
                                                     Column("scaffold_id", Integer, ForeignKey("plant_tribes_scaffold.id"), index=True, nullable=False),
                                                     Column("taxa_id", Integer, ForeignKey("plant_tribes_taxa.id"), index=True, nullable=False),
                                                     Column("orthogroup_id", Integer, ForeignKey("plant_tribes_orthogroup.id"), index=True, nullable=False),
                                                     Column("gene_id", Integer, ForeignKey("plant_tribes_gene.id"), index=True, nullable=False))


mapper(plant_tribes_model.PlantTribesScaffold, plant_tribes_model.PlantTribesScaffold.table)

mapper(plant_tribes_model.PlantTribesTaxa, plant_tribes_model.PlantTribesTaxa.table,
       properties=dict(scaffolds=relation(plant_tribes_model.TaxaScaffoldAssociation,
                                          primaryjoin=(plant_tribes_model.PlantTribesTaxa.table.c.id == plant_tribes_model.TaxaScaffoldAssociation.table.c.taxa_id)),
                       genes=relation(plant_tribes_model.TaxaGeneAssociation,
                                      primaryjoin=(plant_tribes_model.PlantTribesTaxa.table.c.id == plant_tribes_model.TaxaGeneAssociation.table.c.taxa_id))))

mapper(plant_tribes_model.PlantTribesOrthogroup, plant_tribes_model.PlantTribesOrthogroup.table)

mapper(plant_tribes_model.PlantTribesGene, plant_tribes_model.PlantTribesGene.table)

mapper(plant_tribes_model.TaxaScaffoldAssociation, plant_tribes_model.TaxaScaffoldAssociation.table)

mapper(plant_tribes_model.TaxaGeneAssociation, plant_tribes_model.TaxaGeneAssociation.table)

mapper(plant_tribes_model.ScaffoldAssociation, plant_tribes_model.ScaffoldAssociation.table)


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

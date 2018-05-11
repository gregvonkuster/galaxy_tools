import datetime
import logging

from sqlalchemy import Column, DateTime, ForeignKey, Integer, MetaData, Table, TEXT

# Need our custom types, but don't import anything else from model
from galaxy.model.custom_types import TrimmedString

now = datetime.datetime.utcnow
log = logging.getLogger(__name__)
metadata = MetaData()

PlantTribesScaffold_table = Table("plant_tribes_scaffold", metadata,
    Column("scaffold_id", TrimmedString(0), primary_key=True),
    Column("description", TEXT, nullable=False),
    Column("num_genes", Integer, nullable=False),
    Column("num_orthogroups", Integer))


PlantTribesTaxa_table = Table("plant_tribes_taxa", metadata,
    Column("scaffold_id", TrimmedString(0), ForeignKey("plant_tribes_scaffold.scaffold_id"), index=True, nullable=False),
    Column("num_genes", Integer, nullable=False),
    Column("species", TrimmedString(50), index=True, nullable=False),
    Column("species_family", TrimmedString(50), nullable=False),
    Column("species_order", TrimmedString(50), nullable=False),
    Column("species_group", TrimmedString(50), nullable=False),
    Column("species_clade", TrimmedString(50), nullable=False))


PlantTribesOrthogroup_table = Table("plant_tribes_orthogroup", metadata,
    Column("ortho_id", Integer, primary_key=True, nullable=False),
    Column("scaffold_id", TrimmedString(0), ForeignKey("plant_tribes_scaffold.scaffold_id"), index=True, nullable=False),
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

PlantTribesGene_table = Table("plant_tribes_gene", metadata,
    Column("gene_id", Integer, primary_key=True, nullable=False),
    Column("species", TrimmedString(50), nullable=False),
    Column("ortho_id", Integer, ForeignKey("plant_tribes_orthogroup.ortho_id"), index=True, nullable=False),
    Column("scaffold_id", TrimmedString(10), ForeignKey("plant_tribes_scaffold.scaffold_id"), index=True, nullable=False),
    Column("dna_sequence", TEXT, nullable=False),
    Column("aa_sequence", TEXT, nullable=False))

"""
TaxaScaffoldAssociation_table = Table("taxa_scaffold_association", metadata,
    Column("id", Integer, primary_key=True),
    Column("taxa_id", Integer, ForeignKey("plant_tribes_taxa.id"), index=True, nullable=False),
    Column("scaffold_id", Integer, ForeignKey("plant_tribes_scaffold.id"), index=True, nullable=False))

TaxaGeneAssociation_table = Table("taxa_gene_association", metadata,
    Column("id", Integer, primary_key=True),
    Column("taxa_id", Integer, ForeignKey("plant_tribes_taxa.id"), index=True, nullable=False),
    Column("gene_id", Integer, ForeignKey("plant_tribes_gene.id"), index=True, nullable=False))

ScaffoldAssociation_table = Table("scaffold_taxa_orthogroup_gene_association", metadata,
    Column("id", Integer, primary_key=True),
    Column("scaffold_id", Integer, ForeignKey("plant_tribes_scaffold.id"), index=True, nullable=False),
    Column("taxa_id", Integer, ForeignKey("plant_tribes_taxa.id"), index=True, nullable=False),
    Column("orthogroup_id", Integer, ForeignKey("plant_tribes_orthogroup.id"), index=True, nullable=False),
    Column("gene_id", Integer, ForeignKey("plant_tribes_gene.id"), index=True, nullable=False))
"""

def upgrade(migrate_engine):
    print(__doc__)

    metadata.bind = migrate_engine
    metadata.create_all()

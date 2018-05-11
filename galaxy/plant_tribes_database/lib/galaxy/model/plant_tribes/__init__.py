import logging

from galaxy.util.dictifiable import Dictifiable

log = logging.getLogger(__name__)


class PlantTribesScaffold(Dictifiable):
    dict_collection_visible_keys = ['scaffold_id', 'description', 'num_genes', 'num_orthogroups']

    def __init__(self, id=None, description=None, scaffold=None, num_genes=None, num_orthogroups=None):
        self.scaffold_id = scaffold_id
        self.description = description
        self.num_genes = num_genes
        self.num_orthogroups = num_orthogroups

    def as_dict(self, value_mapper=None):
        return self.to_dict(view='element', value_mapper=value_mapper)

    def to_dict(self, view='collection', value_mapper=None):
        if value_mapper is None:
            value_mapper = {}
        rval = {}
        try:
            visible_keys = self.__getattribute__('dict_' + view + '_visible_keys')
        except AttributeError:
            raise Exception('Unknown API view: %s' % view)
        for key in visible_keys:
            try:
                rval[key] = self.__getattribute__(key)
                if key in value_mapper:
                    rval[key] = value_mapper.get(key, rval[key])
            except AttributeError:
                rval[key] = None
        return rval


class PlantTribesTaxa(Dictifiable):
    dict_collection_visible_keys = ['scaffold_id', 'num_genes', 'species', 'species_family', 'species_order', 'species_group', 'species_clade']

    def __init__(self, scaffold_id=None, num_genes=None, species=None, species_family=None, species_order=None, species_group=None, species_clade=None):
        self.scaffold_id = scaffold_id
        self.num_genes = num_genes
        self.species = species
        self.species_family = species_family
        self.species_order = species_order
        self.species_group = species_group
        self.species_clade = species_clade

    def as_dict(self, value_mapper=None):
        return self.to_dict(view='element', value_mapper=value_mapper)

    def to_dict(self, view='collection', value_mapper=None):
        if value_mapper is None:
            value_mapper = {}
        rval = {}
        try:
            visible_keys = self.__getattribute__('dict_' + view + '_visible_keys')
        except AttributeError:
            raise Exception('Unknown API view: %s' % view)
        for key in visible_keys:
            try:
                rval[key] = self.__getattribute__(key)
                if key in value_mapper:
                    rval[key] = value_mapper.get(key, rval[key])
            except AttributeError:
                rval[key] = None
        return rval


class PlantTribesOrthogroup(Dictifiable):
    dict_collection_visible_keys = ['ortho_id', 'scaffold_id', 'num_taxa', 'num_genes', 'super_ortho_1_2', 'super_ortho_1_5', 'super_ortho_1_8',
                                    'super_ortho_2_0', 'super_ortho_2_5', 'super_ortho_3_0', 'super_ortho_3_5', 'super_ortho_4_0',
                                    'super_ortho_4_5', 'super_ortho_5_0', 'ahdr_description', 'tair_description', 'pfam_description',
                                    'interproscan__description', 'gene_ontology_molecular_funcion_description', 'gene_ontology_biological_process_description',
                                    'gene_ontology_celular_component_description']

    def __init__(self, ortho_id=None, scaffold_id=None, num_taxa=None, num_genes=None, super_ortho_1_2=None, super_ortho_1_5=None,
                 super_ortho_1_8=None, super_ortho_2_0=None, super_ortho_2_5=None, super_ortho_3_0=None, super_ortho_3_5=None,
                 super_ortho_4_0=None, super_ortho_4_5=None, super_ortho_5_0=None, ahdr_description=None, tair_description=None,
                 pfam_description=None, interproscan__description=None, gene_ontology_molecular_funcion_description=None,
                 gene_ontology_biological_process_description=None, gene_ontology_celular_component_description=None):
        self.ortho_id = ortho_id
        self.scaffold_id = scaffold_id
        self.num_taxa = num_taxa
        self.num_genes = num_genes
        self.super_ortho_1_2 = super_ortho_1_2
        self.super_ortho_1_5 = super_ortho_1_5
        self.super_ortho_1_8 = super_ortho_1_8
        self.super_ortho_2_0 = super_ortho_2_0
        self.super_ortho_2_5 = super_ortho_2_5
        self.super_ortho_3_0 = super_ortho_3_0
        self.super_ortho_3_5 = super_ortho_3_5
        self.super_ortho_4_0 = super_ortho_4_0
        self.super_ortho_4_5 = super_ortho_4_5
        self.super_ortho_5_0 = super_ortho_5_0
        self.ahdr_description = ahdr_description
        self.tair_description = tair_description
        self.pfam_description = pfam_description
        self.interproscan__description = interproscan__description
        self.gene_ontology_molecular_funcion_description = gene_ontology_molecular_funcion_description
        self.gene_ontology_biological_process_description = gene_ontology_biological_process_description
        self.gene_ontology_celular_component_description = gene_ontology_celular_component_description

    def as_dict(self, value_mapper=None):
        return self.to_dict(view='element', value_mapper=value_mapper)

    def to_dict(self, view='collection', value_mapper=None):
        if value_mapper is None:
            value_mapper = {}
        rval = {}
        try:
            visible_keys = self.__getattribute__('dict_' + view + '_visible_keys')
        except AttributeError:
            raise Exception('Unknown API view: %s' % view)
        for key in visible_keys:
            try:
                rval[key] = self.__getattribute__(key)
                if key in value_mapper:
                    rval[key] = value_mapper.get(key, rval[key])
            except AttributeError:
                rval[key] = None
        return rval


class PlantTribesGene(Dictifiable):
    dict_collection_visible_keys = ['gene_id', 'species', 'ortho_id', 'scaffold_id', 'dna_sequence', 'aa_sequence']

    def __init__(self, gene_id=None, species=None, ortho_id=None, scaffold_id=None, dna_sequence=None, aa_sequence=None):
        self.gene_id = gene_id
        self.species = species
        self.ortho_id = ortho_id
        self.scaffold_id = scaffold_id
        self.dna_sequence = dna_sequence
        self.aa_sequence = aa_sequence

    def as_dict(self, value_mapper=None):
        return self.to_dict(view='element', value_mapper=value_mapper)

    def to_dict(self, view='collection', value_mapper=None):
        if value_mapper is None:
            value_mapper = {}
        rval = {}
        try:
            visible_keys = self.__getattribute__('dict_' + view + '_visible_keys')
        except AttributeError:
            raise Exception('Unknown API view: %s' % view)
        for key in visible_keys:
            try:
                rval[key] = self.__getattribute__(key)
                if key in value_mapper:
                    rval[key] = value_mapper.get(key, rval[key])
            except AttributeError:
                rval[key] = None
        return rval

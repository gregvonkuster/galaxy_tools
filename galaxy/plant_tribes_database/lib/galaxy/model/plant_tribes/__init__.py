import logging

from galaxy.util.dictifiable import Dictifiable

log = logging.getLogger(__name__)


class PlantTribesScaffold(Dictifiable):
    dict_collection_visible_keys = ['id', 'scaffold_id', 'clustering_method']

    def __init__(self, id=None, scaffold_id=None, clustering_method=None):
        self.id = id
        self.scaffold_id = scaffold_id
        self.clustering_method = clustering_method

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


class PlantTribesTaxon(Dictifiable):
    dict_collection_visible_keys = ['id', 'species_name', 'scaffold_id', 'num_genes', 'species_family', 'species_order', 'species_group', 'species_clade']

    def __init__(self, id=None, species_name=None, scaffold_id=None, species_family=None, species_order=None, species_group=None, species_clade=None, num_genes=None):
        self.id = id,
        self.species_name = species_name
        self.scaffold_id = scaffold_id
        self.num_genes = num_genes
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
    dict_collection_visible_keys = ['id', 'orthogroup_id', 'scaffold_id', 'num_species', 'num_genes', 'super_ortho_1_2',
                                    'super_ortho_1_5', 'super_ortho_1_8', 'super_ortho_2_0', 'super_ortho_2_5', 'super_ortho_3_0',
                                    'super_ortho_3_5', 'super_ortho_4_0', 'super_ortho_4_5', 'super_ortho_5_0', 'ahrd_description',
                                    'tair_description', 'pfam_description', 'interproscan_description', 'molecular_function',
                                    'biological_process', 'cellular_component']

    def __init__(self, id=None, orthogroup_id=None, scaffold_id=None, num_species=None, num_genes=None, super_ortho_1_2=None,
                 super_ortho_1_5=None, super_ortho_1_8=None, super_ortho_2_0=None, super_ortho_2_5=None, super_ortho_3_0=None,
                 super_ortho_3_5=None, super_ortho_4_0=None, super_ortho_4_5=None, super_ortho_5_0=None, ahrd_description=None,
                 tair_description=None, pfam_description=None, interproscan_description=None, molecular_function=None, biological_process=None,
                 cellular_component=None):
        self.id = id
        self.orthogroup_id = orthogroup_id
        self.scaffold_id = scaffold_id
        self.num_species = num_species
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
        self.ahrd_description = ahrd_description
        self.tair_description = tair_description
        self.pfam_description = pfam_description
        self.interproscan_description = interproscan_description
        self.molecular_function = molecular_function
        self.biological_process = biological_process
        self.cellular_component = cellular_component

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
    dict_collection_visible_keys = ['id', 'gene_id', 'dna_sequence', 'aa_sequence']

    def __init__(self, id=None, gene_id=None, dna_sequence=None, aa_sequence=None):
        self.id = id,
        self.gene_id = gene_id
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


class GeneScaffoldOrthogroupTaxonAssociation(object):
    dict_collection_visible_keys = ['id', 'gene_id', 'scaffold_id', 'orthogroup_id', 'taxon_id']

    def __init__(self, id=None, gene_id=None, scaffold_id=None, taxon_id=None):
        self.id = id
        self.gene_id = gene_id
        self.scaffold_id = scaffold_id
        self.orthogroup_id = orthogroup_id
        self.taxon_id = taxon_id

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

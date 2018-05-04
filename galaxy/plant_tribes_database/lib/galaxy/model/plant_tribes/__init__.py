import logging

from galaxy.util.dictifiable import Dictifiable

log = logging.getLogger(__name__)


class PlantTribesScaffold(Dictifiable):
    dict_collection_visible_keys = ['id', 'description', 'genes', 'orthogroups']

    def __init__(self, id=None, create_time=None, description=None, genes=None, orthogroups=None):
        self.id = id
        self.create_time = create_time
        self.description = description
        self.genes = genes
        self.orthogroups = orthogroups

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
    dict_collection_visible_keys = ['id', 'genes', 'species_family', 'species_order', 'species_group', 'species_clade']

    def __init__(self, id=None, create_time=None, genes=None, species_family=None, species_order=None, species_group=None, species_clade=None):
        self.id = id
        self.create_time = create_time
        self.genes = genes
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
    dict_collection_visible_keys = ['id', 'taxa', 'genes', 'super_ortho_1_2', 'super_ortho_1_5', 'super_ortho_1_8', 'super_ortho_2_0', 'super_ortho_2_5',
                                    'super_ortho_3_0', 'super_ortho_3_5', 'super_ortho_4_0', 'super_ortho_4_5', 'super_ortho_5_0', 'ahdr_desc', 'tair_desc',
                                    'pfam_desc', 'iprscan_desc', 'go_mf_desc', 'go_bp_desc', 'go_cc_desc']

    def __init__(self, id=None, create_time=None, taxa=None, genes=None, super_ortho_1_2=None, super_ortho_1_5=None, super_ortho_1_8=None, super_ortho_2_0=None,
                 super_ortho_2_5=None, super_ortho_3_0=None, super_ortho_3_5=None, super_ortho_4_0=None, super_ortho_4_5=None, super_ortho_5_0=None,
                 ahdr_desc=None, tair_desc=None, pfam_desc=None, iprscan_desc=None, go_mf_desc=None, go_bp_desc=None, go_cc_desc=None):
        self.id = id
        self.create_time = create_time
        self.taxa = taxa
        self.genes = genes
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
        self.ahdr_desc = ahdr_desc
        self.tair_desc = tair_desc
        self.pfam_desc = pfam_desc
        self.iprscan_desc = iprscan_desc
        self.go_mf_desc = go_mf_desc
        self.go_bp_desc = go_bp_desc
        self.go_cc_desc = go_cc_desc

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
    dict_collection_visible_keys = ['id', 'dna_sequence', 'aa_sequence']

    def __init__(self, id=None, create_time=None, dna_sequence=None, aa_sequence=None):
        self.id = id
        self.create_time = create_time
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


class ScaffoldAssociation(Dictifiable):
    dict_collection_visible_keys = ['id', 'scaffold', 'scaffold_id', 'taxa_id', 'orthogroup_id', 'gene_id']

    def __init__(self, id=None, create_time=None, scaffold=None, scaffold_id=None, taxa_id=None, orthogroup_id=None, gene_id=None):
        self.id = id
        self.create_time = create_time
        self.scaffold = scaffold
        self.scaffold_id = scaffold_id
        self.taxa_id = taxa_id
        self.orthogroup_id = orthogroup_id
        self.gene_id = gene_id


class SpeciesAssociation(Dictifiable):
    dict_collection_visible_keys = ['id', 'species', 'taxa_id', 'gene_id']

    def __init__(self, id=None, create_time=None, species=None, taxa_id=None, gene_id=None):
        self.id = id
        self.create_time = create_time
        self.species = species
        self.taxa_id = taxa_id
        self.gene_id = gene_id

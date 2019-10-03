import logging

from galaxy.util.dictifiable import Dictifiable

log = logging.getLogger(__name__)


class Allele(Dictifiable):
    dict_collection_visible_keys = ['id', 'allele']

    def __init__(self, allele=None):
        self.allele = allele

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


class Collector(Dictifiable):
    dict_collection_visible_keys = ['id', 'person_id', 'contact_id']

    def __init__(self, person_id=None, contact_id=None):
        self.person_id = person_id
        self.contact_id = contact_id

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


class Colony(Dictifiable):
    dict_collection_visible_keys = ['id', 'latitude', 'longitude', 'depth', 'reef_id']

    def __init__(self, latitude=None, longitude=None, depth=None, reef_id=None):
        # Latitude in decimal degrees, WGS84 datum, geographic location of
        # the wild donor colony from which the sample was collected.
        self.latitude = latitude
        # Longitude in decimal degrees, WGS84 datum, geographic location of
        # the wild donor colony from which the sample was collected.
        self.longitude = longitude
        # Depth in meters of the wild donor colony.
        self.depth = depth
        self.reef_id = reef_id

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


class Experiment(Dictifiable):
    dict_collection_visible_keys = ['id', 'seq_facility', 'array_version']

    def __init__(self, seq_facility=None, array_version=None):
        # Description of experiment metadata.
        self.seq_facility = seq_facility
        self.array_version = array_version

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


class Fragment(Dictifiable):
    dict_collection_visible_keys = ['id', 'colony_id']

    def __init__(self, colony_id=None):
        self.colony_id = colony_id

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


class Genotype(Dictifiable):
    dict_collection_visible_keys = ['id', 'coral_mlg_clonal_id', 'coral_mlg_rep_sample_id', 'genetic_coral_species_call']

    def __init__(self, coral_mlg_clonal_id=None, coral_mlg_rep_sample_id=None, genetic_coral_species_call=None):
        # Description of experiment metadata.
        self.coral_mlg_clonal_id = coral_mlg_clonal_id
        self.coral_mlg_rep_sample_id = coral_mlg_rep_sample_id
        self.genetic_coral_species_call = genetic_coral_species_call

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


class Person(Dictifiable):
    dict_collection_visible_keys = ['id', 'last_name', 'first_name', 'organization', 'email']

    def __init__(self, last_name=None, first_name=None, organization=None, email=None):
        self.last_name = last_name
        self.first_name = first_name
        self.organization = organization
        self.email = email

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


class ProbeAnnotation(Dictifiable):
    dict_collection_visible_keys = ['id', 'probe_set_id', 'affy_snp_id', 'chr_id', 'start',
                                    'strand', 'flank', 'allele_a', 'allele_b', 'allele_frequencies',
                                    'annotation_notes', 'allele_count', 'ordered_alleles', 'chrtype', 'custchr',
                                    'custid', 'custpos', 'organism', 'pconvert', 'recommendation',
                                    'refstr', 'snppriority', 'genotype_probe', 'fixed_status', 'acerv_allele']

    def __init__(self, probe_set_id=None, affy_snp_id=None, chr_id=None, start=None,
                 strand=None, flank=None, allele_a=None, allele_b=None, allele_frequencies=None,
                 annotation_notes=None, allele_count=None, ordered_alleles=None, chrtype=None, custchr=None,
                 custid=None, custpos=None, organism=None, pconvert=None, recommendation=None,
                 refstr=None, snppriority=None, genotype_probe=None, fixed_status=None, acerv_allele=None):
        self.probe_set_id = probe_set_id
        self.affy_snp_id = affy_snp_id
        self.chr_id = chr_id
        self.start = start
        self.strand = strand
        self.flank = flank
        self.allele_a = allele_a
        self.allele_b = allele_b
        self.allele_frequencies = allele_frequencies
        self.annotation_notes = annotation_notes
        self.allele_count = allele_count
        self.ordered_alleles = ordered_alleles
        self.chrtype = chrtype
        self.custchr = custchr
        self.custid = custid
        self.custpos = custpos
        self.organism = organism
        self.pconvert = pconvert
        self.recommendation = recommendation
        self.refstr = refstr
        self.snppriority = snppriority
        self.genotype_probe = genotype_probe
        self.fixed_status = fixed_status
        self.acerv_allele = acerv_allele

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


class Reef(Dictifiable):
    dict_collection_visible_keys = ['id', 'name', 'region', 'latitude', 'longitude', 'geographic_origin']

    def __init__(self, name=None, region=None, latitude=None, longitude=None, geographic_origin=None):
        self.name = name
        self.region = region
        # Latitude in decimal degrees, WGS84 datum, geographic location of reef.
        self.latitude = latitude
        # Longitude in decimal degrees, WGS84 datum, geographic location of reef.
        self.longitude = longitude
        self.geographic_origin = geographic_origin

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


class Phenotype(Dictifiable):
    dict_collection_visible_keys = ['id', 'disease_resist', 'bleach_resist', 'mortality', 'tle',
                                    'spawning', 'sperm_motility', 'healing_time']

    def __init__(self, disease_resist=None, bleach_resist=None, mortality=None, tle=None,
                 spawning=None, sperm_motility=None, healing_time=None):
        # Description of experiment metadata.
        self.disease_resist = disease_resist
        self.bleach_resist = bleach_resist
        self.mortality = mortality
        self.tle = tle
        self.spawning = spawning
        self.sperm_motility = sperm_motility
        self.healing_time = healing_time

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


class Sample(Dictifiable):
    dict_collection_visible_keys = ['id', 'create_time', 'affy_id', 'sample_id', 'allele_id',
                                    'genotype_id', 'phenotype_id', 'experiment_id', 'colony_id', 'colony_location',
                                    'fragment_id', 'taxonomy_id', 'collector_id', 'collection_date', 'user_specimen_id',
                                    'registry_id', 'depth', 'dna_extraction_method', 'dna_concentration', 'public',
                                    'public_after_date,' 'percent_missing_data_coral', 'percent_missing_data_sym',
                                    'percent_acerv_coral', 'percent_reference_sym', 'percent_apalm_coral',
                                    'percent_alternative_sym', 'percent_heterozygous_coral', 'percent_heterozygous_sym',
                                    'field_call', 'bcoral_genet_id']

    def __init__(self, create_time=None, affy_id=None, sample_id=None, allele_id=None,
                 genotype_id=None, phenotype_id=None, experiment_id=None, colony_id=None, colony_location=None,
                 fragment_id=None, taxonomy_id=None, collector_id=None, collection_date=None, user_specimen_id=None,
                 registry_id=None, depth=None, dna_extraction_method=None, dna_concentration=None, public=None,
                 public_after_date=None, percent_missing_data_coral=None, percent_missing_data_sym=None,
                 percent_acerv_coral=None, percent_reference_sym=None, percent_apalm_coral=None,
                 percent_alternative_sym=None, percent_heterozygous_coral=None, percent_heterozygous_sym=None,
                 field_call=None, bcoral_genet_id=None):
        self.create_time = create_time
        self.affy_id = affy_id
        self.sample_id = sample_id
        self.allele_id = allele_id
        self.genotype_id = genotype_id
        self.phenotype_id = phenotype_id
        self.experiment_id = experiment_id
        self.colony_id = colony_id
        self.colony_location = colony_location
        self.fragment_id = fragment_id
        self.taxonomy_id = taxonomy_id
        self.collector_id = collector_id
        self.collection_date = collection_date
        # User specific identification for samples, should limit to some sort of length.
        self.user_specimen_id = user_specimen_id
        self.registry_id = registry_id
        # Depth in meters of the sample for this entry.
        self.depth = depth
        # Method used to extract DNA.
        self.dna_extraction_method = dna_extraction_method
        # DNA concentration in ug.
        self.dna_concentration = dna_concentration
        self.public = public
        self.public_after_date = public_after_date
        self.percent_missing_data_coral = percent_missing_data_coral
        self.percent_missing_data_sym = percent_missing_data_sym
        self.percent_acerv_coral = percent_acerv_coral
        self.percent_reference_sym = percent_reference_sym
        self.percent_apalm_coral = percent_apalm_coral
        self.percent_alternative_sym = percent_alternative_sym
        self.percent_heterozygous_coral = percent_heterozygous_coral
        self.percent_heterozygous_sym = percent_heterozygous_sym
        self.field_call = field_call
        self.bcoral_genet_id = bcoral_genet_id

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


class Taxonomy(Dictifiable):
    dict_collection_visible_keys = ['id', 'species_name', 'genus_name']

    def __init__(self, species_name=None, genus_name=None):
        self.species_name = species_name
        self.genus_name = genus_name

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

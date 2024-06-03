import datetime
import logging
from typing import TYPE_CHECKING

from sqlalchemy import (
    Boolean,
    Column,
    DateTime,
    ForeignKey,
    Integer,
    Numeric,
    String,
    Text,
)
from sqlalchemy.orm import (
    registry,
    relationship,
)
from galaxy.model.custom_types import TrimmedString
from galaxy.model.orm.now import now
from galaxy.util.dictifiable import Dictifiable

if TYPE_CHECKING:
    # Workaround for https://github.com/python/mypy/issues/14182
    from sqlalchemy.orm.decl_api import DeclarativeMeta as _DeclarativeMeta

    class DeclarativeMeta(_DeclarativeMeta, type):
        pass

else:
    from sqlalchemy.orm.decl_api import DeclarativeMeta

log = logging.getLogger(__name__)

mapper_registry = registry()

# Get current date plus one year for insertion into
# the public_after_date column of the sample table.
today = datetime.date.today()
try:
    # Return the same day of the year.
    year = today.year + 1
    year_from_now = today.replace(year=year)
except Exception:
    # Handle leap years.
    year_from_now = today + (datetime.date(today.year + 1, 1, 1) - datetime.date(today.year, 1, 1))


def cached_id(corals_model_object):
    """Get model object id attribute without a firing a database query.

    Useful to fetching the id of a typical corals model object after a flush,
    where SA is going to mark the id attribute as unloaded but we know the id
    is immutable and so we can use the database identity to fetch.

    With Galaxy's default SA initialization - any flush marks all attributes as
    unloaded - even objects completely unrelated to the flushed changes and
    even attributes we know to be immutable like id. See test_galaxy_mapping.py
    for verification of this behavior. This method is a workaround that uses
    the fact that we know all corals objects use the id attribute as identity
    and SA internals (_sa_instance_state) to infer the previously loaded ID
    value. I tried digging into the SA internals extensively and couldn't find
    a way to get the previously loaded values after a flush to allow a
    generalization of this for other attributes.
    """
    if hasattr(corals_model_object, "_sa_instance_state"):
        identity = corals_model_object._sa_instance_state.identity
        if identity:
            assert len(identity) == 1
            return identity[0]

    return corals_model_object.id


class Base(metaclass=DeclarativeMeta):
    __abstract__ = True
    registry = mapper_registry
    metadata = mapper_registry.metadata
    __init__ = mapper_registry.constructor

    @classmethod
    def __declare_last__(cls):
        cls.table = cls.__table__


class RepresentById:
    id: int

    def __repr__(self):
        try:
            r = f"<galaxy.model.corals.{self.__class__.__name__}({cached_id(self)}) at {hex(id(self))}>"
        except Exception:
            r = object.__repr__(self)
            log.exception("Caught exception attempting to generate repr for: %s", r)
        return r


class Allele(Base, Dictifiable, RepresentById):
    __tablename__ = "allele"

    id = Column(Integer, primary_key=True)
    create_time = Column(DateTime, default=now)
    update_time = Column(DateTime, default=now, onupdate=now)
    allele = Column("allele", Text)
    samples = relationship("Sample", back_populates="alleles")

    dict_element_visible_keys = ['id', 'allele']

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


class Colony(Base, Dictifiable, RepresentById):
    __tablename__ = "colony"

    id = Column(Integer, primary_key=True)
    create_time = Column(DateTime, default=now)
    update_time = Column(DateTime, default=now, onupdate=now)
    latitude = Column(Numeric(15, 6))
    longitude = Column(Numeric(15, 6))
    depth = Column(Numeric(15, 6))
    reef_id = Column(Integer, ForeignKey("reef.id"), index=True)
    reef = relationship("Reef", back_populates="colonies")
    fragments = relationship("Fragment", back_populates="colony")
    matching_samples = relationship("Sample", back_populates="colony")

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


class Experiment(Base, Dictifiable, RepresentById):
    __tablename__ = "experiment"

    id = Column(Integer, primary_key=True)
    create_time = Column(DateTime, default=now)
    update_time = Column(DateTime, default=now, onupdate=now)
    seq_facility = Column(String)
    array_version = Column(TrimmedString(255))
    result_folder_name = Column(TrimmedString(255))
    plate_barcode = Column(TrimmedString(255))
    samples = relationship("Sample", back_populates="experiment")

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


class Fragment(Base, Dictifiable, RepresentById):
    __tablename__ = "fragment"

    id = Column(Integer, primary_key=True)
    create_time = Column(DateTime, default=now)
    update_time = Column(DateTime, default=now, onupdate=now)
    colony_id = Column(Integer, ForeignKey("colony.id"), index=True)
    colony = relationship("Colony", back_populates="fragments")

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


class Genotype(Base, Dictifiable, RepresentById):
    __tablename__ = "genotype"

    id = Column(Integer, primary_key=True)
    create_time = Column(DateTime, default=now)
    update_time = Column(DateTime, default=now, onupdate=now)
    coral_mlg_clonal_id = Column(TrimmedString(255))
    coral_mlg_rep_sample_id = Column(TrimmedString(255))
    genetic_coral_species_call = Column(TrimmedString(255))
    samples = relationship("Sample", back_populates="genotype")

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


class Person(Base, Dictifiable, RepresentById):
    __tablename__ = "person"

    id = Column(Integer, primary_key=True)
    create_time = Column(DateTime, default=now)
    update_time = Column(DateTime, default=now, onupdate=now)
    last_name = Column(TrimmedString(255))
    first_name = Column(TrimmedString(255))
    organization = Column(TrimmedString(255))
    email = Column(TrimmedString(255))
    samples = relationship("Sample", back_populates="collector")

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


class ProbeAnnotation(Base, Dictifiable, RepresentById):
    __tablename__ = "probe_annotation"

    id = Column(Integer, primary_key=True)
    create_time = Column(DateTime, default=now)
    update_time = Column(DateTime, default=now, onupdate=now)
    probe_set_id = Column(TrimmedString(255))
    affy_snp_id = Column(TrimmedString(255))
    chr_id = Column(Integer)
    start = Column(Integer)
    strand = Column(TrimmedString(255))
    flank = Column(TrimmedString(255))
    allele_a = Column(TrimmedString(255))
    allele_b = Column(TrimmedString(255))
    allele_frequencies = Column(TrimmedString(255))
    annotation_notes = Column(TrimmedString(255))
    allele_count = Column(TrimmedString(255))
    ordered_alleles = Column(TrimmedString(255))
    chrtype = Column(TrimmedString(255))
    custchr = Column(TrimmedString(255))
    custid = Column(TrimmedString(255))
    custpos = Column(TrimmedString(255))
    organism = Column(TrimmedString(255))
    pconvert = Column(TrimmedString(255))
    recommendation = Column(TrimmedString(255))
    refstr = Column(TrimmedString(255))
    snppriority = Column(TrimmedString(255))
    genotype_probe = Column(TrimmedString(255))
    fixed_status = Column(TrimmedString(255))
    acerv_allele = Column(TrimmedString(255))

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


class Reef(Base, Dictifiable, RepresentById):
    __tablename__ = "reef"

    id = Column(Integer, primary_key=True)
    create_time = Column(DateTime, default=now)
    update_time = Column(DateTime, default=now, onupdate=now)
    name = Column(TrimmedString(255))
    region = Column(TrimmedString(255))
    latitude = Column(Numeric(15, 6))
    longitude = Column(Numeric(15, 6))
    geographic_origin = Column(TrimmedString(255))
    colonies = relationship("Colony", back_populates="reef")

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


class Phenotype(Base, Dictifiable, RepresentById):
    __tablename__ = "phenotype"

    id = Column(Integer, primary_key=True)
    create_time = Column(DateTime, default=now)
    update_time = Column(DateTime, default=now, onupdate=now)
    disease_resist = Column(TrimmedString(255))
    bleach_resist = Column(TrimmedString(255))
    mortality = Column(TrimmedString(255))
    tle = Column(TrimmedString(255))
    spawning = Column(Boolean)
    sperm_motility = Column(Numeric(15, 6))
    healing_time = Column(Numeric(15, 6))
    samples = relationship("Sample", back_populates="phenotype")

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


class Sample(Base, Dictifiable, RepresentById):
    __tablename__ = "sample"

    id = Column(Integer, primary_key=True)
    create_time = Column(DateTime, default=now)
    update_time = Column(DateTime, default=now, onupdate=now)
    affy_id = Column(TrimmedString(255), index=True, nullable=False)
    sample_id = Column(TrimmedString(255), index=True, nullable=False)
    allele_id = Column(Integer, ForeignKey("allele.id"), index=True)
    genotype_id = Column(Integer, ForeignKey("genotype.id"), index=True)
    phenotype_id = Column(Integer, ForeignKey("phenotype.id"), index=True)
    experiment_id = Column(Integer, ForeignKey("experiment.id"), index=True)
    colony_id = Column(Integer, ForeignKey("colony.id"), index=True)
    colony_location = Column(TrimmedString(255))
    taxonomy_id = Column(Integer, ForeignKey("taxonomy.id"), index=True)
    collector_id = Column(Integer, ForeignKey("person.id"), index=True)
    collection_date = Column(DateTime)
    user_specimen_id = Column(TrimmedString(255))
    registry_id = Column(TrimmedString(255))
    depth = Column(Numeric(15, 6))
    dna_extraction_method = Column(TrimmedString(255))
    dna_concentration = Column(Numeric(10, 6))
    public = Column(Boolean)
    public_after_date = Column(DateTime, default=year_from_now)
    percent_missing_data_coral = Column(Numeric(15, 6))
    percent_missing_data_sym = Column(Numeric(15, 6))
    percent_acerv_coral = Column(Numeric(15, 6))
    percent_reference_sym = Column(Numeric(15, 6))
    percent_apalm_coral = Column(Numeric(15, 6))
    percent_alternative_sym = Column(Numeric(15, 6))
    percent_heterozygous_coral = Column(Numeric(15, 6))
    percent_heterozygous_sym = Column(Numeric(15, 6))
    field_call = Column(TrimmedString(255))
    bcoral_genet_id = Column(TrimmedString(255))
    alleles = relationship("Allele", back_populates="samples")
    genotype = relationship("Genotype", back_populates="samples")
    phenotype = relationship("Phenotype", back_populates="samples")
    experiment = relationship("Experiment", back_populates="samples")
    colony = relationship("Colony", back_populates="matching_samples")
    taxonomy = relationship("Taxonomy", back_populates="samples")
    collector = relationship("Person", back_populates="samples")

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


class Taxonomy(Base, Dictifiable, RepresentById):
    __tablename__ = "taxonomy"

    id = Column(Integer, primary_key=True)
    create_time = Column(DateTime, default=now)
    update_time = Column(DateTime, default=now, onupdate=now)
    species_name = Column(TrimmedString(255))
    genus_name = Column(TrimmedString(255))
    samples = relationship("Sample", back_populates="taxonomy")

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

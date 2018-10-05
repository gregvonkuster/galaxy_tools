
from datetime import datetime
from sqlalchemy import (
    and_,
    asc,
    Boolean,
    Column,
    DateTime,
    desc,
    false,
    ForeignKey,
    func,
    Integer,
    MetaData,
    not_,
    Numeric,
    select,
    String,
    Table,
    TEXT,
    Text,
    true,
    Unicode,
    UniqueConstraint,
    VARCHAR
)
from sqlalchemy.ext.associationproxy import association_proxy
from sqlalchemy.ext.orderinglist import ordering_list
from sqlalchemy.orm import backref, class_mapper, column_property, deferred, mapper, object_session, relation
from sqlalchemy.orm.collections import attribute_mapped_collection
from sqlalchemy.sql import exists
from sqlalchemy.types import BigInteger
from sqlalchemy import MetaData
from sqlalchemy import create_engine

now = datetime.utcnow
database_file = "coralSNPs.sqlite"
database_connection = "sqlite:///./%s?isolation_level=IMMEDIATE" % database_file
engine = create_engine(database_connection)
metadata = MetaData(bind=engine)

class Coralvcf_alleles():

    def __init__(self, sample_id=None):
        self.sample_id = sample_id


Coralvcf_alleles.table = Table("coralvcf_alleles", metadata,
    Column("id", Integer, primary_key=True),
    Column("create_time", DateTime, default=now),
    Column("update_time", DateTime, default=now, onupdate=now),
    Column("sample_id", String, index=True, nullable=False))

metadata.create_all()
# conn = engine.connect()
# do more stuff here...


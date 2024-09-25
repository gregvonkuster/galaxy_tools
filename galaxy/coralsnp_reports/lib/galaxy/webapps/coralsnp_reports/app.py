import logging
import sys
import time

import galaxy.model.corals.mapping
from galaxy.config import configure_logging
from galaxy.model.base import SharedModelMapping
from galaxy.security import idencoding
from galaxy.structured_app import BasicSharedApp
from galaxy.web_stack import application_stack_instance
from . import config

log = logging.getLogger(__name__)


class UniverseApplication(BasicSharedApp):
    """Encapsulates the state of a Universe application"""

    def __init__(self, **kwargs):
        super().__init__()
        self[BasicSharedApp] = self
        self.name = "coralsnp_reports"
        # Read config file and check for errors
        self.config = config.Configuration(**kwargs)
        self.config.check()
        configure_logging(self.config)
        log.debug("XXX kwargs: %s\n" % str(kwargs))
        log.debug("XXX self.config: %s\n" % str(self.config))
        log.debug("XXX self.config.corals_database_connection: %s\n" % self.config.corals_database_connection)
        log.debug("python path is: %s", ", ".join(sys.path))
        self.application_stack = application_stack_instance()
        # Determine the database url
        if self.config.corals_database_connection:
            db_url = self.config.corals_database_connection
        else:
            db_url = f"sqlite:///{self.config.database}?isolation_level=IMMEDIATE"
        # Setup the database engine and ORM
        self.model = galaxy.model.corals.mapping.init(self.config.file_path, db_url, self.config.database_engine_options)
        # Security helper
        self.security = idencoding.IdEncodingHelper(id_secret=self.config.id_secret)

        self._register_singleton(idencoding.IdEncodingHelper, self.security)
        self._register_singleton(SharedModelMapping, self.model)

        # used for cachebusting -- refactor this into a *SINGLE* UniverseApplication base.
        self.server_starttime = int(time.time())

    def shutdown(self):
        pass

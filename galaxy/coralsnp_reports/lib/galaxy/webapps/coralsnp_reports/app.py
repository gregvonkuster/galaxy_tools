import logging
import sys
import time

import galaxy.model.corals
import galaxy.model.corals.mapping
from galaxy.security import idencoding
from galaxy.web.stack import application_stack_instance
from . import config

log = logging.getLogger(__name__)


class UniverseApplication(object):
    """Encapsulates the state of a Universe application"""

    def __init__(self, **kwargs):
        log.debug("python path is: %s", ", ".join(sys.path))
        self.name = "coralsnp_reports"
        # Read config file and check for errors
        self.config = config.Configuration(**kwargs)
        self.config.check()
        config.configure_logging(self.config)
        self.application_stack = application_stack_instance()
        # Determine the database url
        if self.config.database_connection:
            db_url = self.config.database_connection
        else:
            db_url = "postgresql://galaxy:galaxy@localhost/stag"
        # Setup the database engine and ORM
        self.model = galaxy.model.corals.mapping.init(db_url,
                                                      self.config.database_engine_options)
        if not self.config.database_connection:
            self.targets_mysql = False
        else:
            self.targets_mysql = 'mysql' in self.config.database_connection
        # Security helper
        self.security = idencoding.IdEncodingHelper(id_secret=self.config.id_secret)
        # used for cachebusting -- refactor this into a *SINGLE* UniverseApplication base.
        self.server_starttime = int(time.time())

    def shutdown(self):
        pass

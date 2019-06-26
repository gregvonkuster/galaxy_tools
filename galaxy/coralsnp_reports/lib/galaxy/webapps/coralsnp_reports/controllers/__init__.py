"""Galaxy coralsnp reports controllers."""
import logging

from galaxy import util

log = logging.getLogger(__name__)


class BaseController(object):
    """
    Base class for CoralSNP Reports web application controllers.
    """

    def __init__(self, app):
        """Initialize an interface for application 'app'"""
        self.app = app
        self.sa_session = app.model.context


class BaseUIController(BaseController):

    def get_object(self, trans, id, class_name, check_ownership=False, check_accessible=False, deleted=None):
        try:
            return BaseController.get_object(self, trans, id, class_name, check_ownership=check_ownership,
                                             check_accessible=check_accessible, deleted=deleted)
        except Exception:
            log.exception("Exception in get_object check for %s %s:", class_name, str(id))
            raise Exception('Server error retrieving %s id ( %s ).' % (class_name, str(id)))

    def message_exception(self, trans, message, sanitize=True):
        trans.response.status = 400
        return {'err_msg': util.sanitize_text(message) if sanitize else message}

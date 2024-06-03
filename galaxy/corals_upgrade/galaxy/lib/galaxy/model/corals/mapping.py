import logging
from threading import local
from typing import (
    Optional,
    Type,
)

from galaxy.model import corals as corals_model
from galaxy.model.base import ModelMapping
from galaxy.model.orm.engine_factory import build_engine
from galaxy.model.corals import mapper_registry

log = logging.getLogger(__name__)

metadata = mapper_registry.metadata


class GalaxyModelMapping(ModelMapping):
    security_agent: None
    thread_local_log: Optional[local]
    GalaxySession: Type


def init(
    file_path,
    url,
    engine_options=None,
    trace_logger=None,
    slow_query_log_threshold=0,
    thread_local_log: Optional[local] = None,
    log_query_counts=False,
) -> GalaxyModelMapping:
    # Build engine
    engine = build_engine(
        url,
        engine_options,
        trace_logger,
        slow_query_log_threshold,
        thread_local_log=thread_local_log,
        log_query_counts=log_query_counts,
    )

    # Configure model, build ModelMapping
    return configure_model_mapping(file_path, engine, thread_local_log)


def configure_model_mapping(
    file_path: str,
    engine,
    thread_local_log,
) -> GalaxyModelMapping:
    return _build_model_mapping(engine, thread_local_log)


def _build_model_mapping(engine, thread_local_log) -> GalaxyModelMapping:
    model_modules = [corals_model]
    model_mapping = GalaxyModelMapping(model_modules, engine)
    model_mapping.thread_local_log = thread_local_log
    return model_mapping
